# import pandas as pd
# import numpy as np
# import statsmodels.formula.api as smf
# from statsmodels.stats.multitest import multipletests, local_fdr
# from scipy.stats import zscore, uniform
# import matplotlib.pyplot as plt
from scipy import stats
import math
import warnings
import requests
import pandas as pd
import statistics
warnings.filterwarnings("ignore")

def mtsFinderEngine(peptides, conditions, pairs, normalization):

    if normalization: #Abundances Normalized F1 126 Sample Baseline #Abundance F1 126 Sample Baseline
        channels = [col for col in peptides.columns if 'Abundances (Normalized)' in col]
        if channels == []:
            channels = [col for col in peptides.columns if 'Abundances Normalized' in col]
    else:
        channels = [col for col in peptides.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in peptides.columns if 'Abundance' in col]

    s = 0
    for condition in conditions: # to remove abundances to skip (empty channels)
        if str(condition.lower()) == 'skip':
            peptides.drop(channels[s], axis=1, inplace=True)
        # elif str(condition.lower()) == 'light':
        #     peptides.drop(channels[s], axis=1, inplace=True)
        elif str(condition.lower()) == 'boost':
            peptides.drop(channels[s], axis=1, inplace=True)
        s+=1

    conditions = [x for x in conditions if not str(x).lower().__contains__('skip')
                  # and not str(x).lower().__contains__('light')
                  and not str(x).lower().__contains__('boost')] # to remove 'skip' from the conditions

    conditions = [i.strip() for i in conditions]
    conditions = [i.lstrip() for i in conditions]

    # to use it second time, because we are dropping columns and second time we are getting column names again after drop
    if normalization: #Abundances Normalized F1 126 Sample Baseline #Abundance F1 126 Sample Baseline
        channels = [col for col in peptides.columns if 'Abundances (Normalized)' in col]
        if channels == []:
            channels = [col for col in peptides.columns if 'Abundances Normalized' in col]
    else:
        channels = [col for col in peptides.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in peptides.columns if 'Abundance' in col]

    columnDict = {channels[i]: conditions[i] for i in range(len(channels))}
    print(columnDict)
    result = MTS_finder(peptides, channels, conditions, pairs)

    result = result.rename(columns=columnDict)

    if pairs[0] != ['']:
        result = calculations(result, conditions, pairs)

    return result

def MTS_finder(data, channels, conditions, pairs):

    columns = ["Accession", "Gene Symbol", "Positions in Master Proteins", "Modifications", "Uniprot", "Uniprot Location",
               "TargetP", "TargetP Location"]

    for channel in channels:  # channels name inserted in columns
        columns.append(channel)

    MTS_uniprot = pd.read_excel('./files/Uniprot_MTS.xlsx')
    MTS_targetP = pd.read_excel('./files/TargetP2_0_prediction_mitocarta3.xlsx')
    MTS_combinedAccession = pd.read_excel('./files/Total_accession_precurser.xlsx')

    accessions_position = list(data['Positions in Master Proteins'])

    out_df = pd.DataFrame(columns=columns)

    num = 0
    for accesion in accessions_position:

        if str(accesion) == 'nan':
            num += 1
            continue

        sum = 0 # finding empty values
        for i in range(0, len(channels)):
            sum = sum + float(data[channels[i]][num])

        if str(sum) == 'nan': # removing empty values
            num += 1
            continue

        if ';' in accesion:
            accesionList = accesion.split('; ')  # double accession ids
            for newAccession in accesionList:
                try:
                    accesionFinal = newAccession.strip().split(' [')[0].split('-')[0]
                except:
                    accesionFinal = newAccession.strip().split(' [')[0]

                if accesionFinal in list(MTS_combinedAccession['Accession']):
                    # first check accession is in the MTS precurser list in the begining before sendin to checker
                    modification = data["Modifications"][num]
                    cache = uniprot(accesion, accesionFinal, modification, MTS_uniprot)

                    if cache != None:
                        for i in range(0, len(channels)):
                            cache[f'{columns[i + 8]}'] = data[channels[i]][num]

                    # out_df = out_df.append(cache, ignore_index=True)
                    out_df = pd.concat([out_df, pd.DataFrame([cache])], ignore_index=True)

                    if cache != None:
                        for i in range(0, len(channels)):
                            cache[f'{columns[i + 8]}'] = data[channels[i]][num]

                    cache = targetP(accesion, accesionFinal, modification, MTS_targetP)
                    # out_df = out_df.append(cache, ignore_index=True)
                    out_df = pd.concat([out_df, pd.DataFrame([cache])], ignore_index=True)

        try:
            accesionFinal = accesion.split(' [')[0].split('-')[0]
        except:
            accesionFinal = accesion.split(' [')[0]

        if accesionFinal in list(MTS_combinedAccession['Accession']):
            modification = data["Modifications"][num]
            cache = uniprot(accesion, accesionFinal, modification, MTS_uniprot)

            if cache != None:
                for i in range(0, len(channels)):
                    cache[f'{columns[i+8]}'] = data[channels[i]][num]

            # out_df = out_df.append(cache, ignore_index=True)
            out_df = pd.concat([out_df, pd.DataFrame([cache])], ignore_index=True)
            cache = targetP(accesion, accesionFinal, modification, MTS_targetP)

            if cache != None:
                for i in range(0, len(channels)):
                    cache[f'{columns[i+8]}'] = data[channels[i]][num]

            # out_df = out_df.append(cache, ignore_index=True)
            out_df = pd.concat([out_df, pd.DataFrame([cache])], ignore_index=True)

        num += 1
    return out_df

def uniprot(accesionwLocation, accesionFinal, modification, MTS_uniprot):

    # check wheter MTS or not!
    entrylist_uniprot = list(MTS_uniprot['Uniprot ID'])
    entryNumberList_uniprot = list(MTS_uniprot['Presequence Location'])
    # P49189 [275-293]
    # Q00610 [912-923]; P53675 [912-923]

    j = 0

    for entry in entrylist_uniprot:
        if accesionFinal == entry:
            break
        else:
            j += 1

    entryNumber = int(entryNumberList_uniprot[j-1])

    # locationLast = int(accesionwLocation.split('[')[1].split('-')[1].split(']')[0])
    locationFirst = int(accesionwLocation.split('[')[1].split('-')[0])

    if locationFirst <= entryNumber:
        try:
            url = f'https://www.ebi.ac.uk/proteins/api/proteins/{accesionFinal}'
            req = requests.get(url)
            result = req.json()
            geneSymbol = result['gene'][0]['name']['value']
        except:
            try:
                url = f'https://www.ebi.ac.uk/proteins/api/proteins/{accesionFinal}'
                req = requests.get(url)
                result = req.json()
                geneSymbol = result['gene'][0]['name']['value']
            except:
                geneSymbol = accesionFinal

        geneSymbolModification = f"{geneSymbol} [{accesionwLocation.split('[')[1]}"

        try:
            geneSymbolModification = geneSymbolModification.split(';')[0]
        except:
            pass

        uniprotDecision = 'Yes'

        cache = {
            "Accession": accesionFinal,
            "Gene Symbol": geneSymbol,
            "Modifications": modification,
            "Positions in Master Proteins": geneSymbolModification,
            "Uniprot": uniprotDecision,
            "Uniprot Location": entryNumber,
            "TargetP": 'No',
            "TargetP Location": 0,
        }

        return cache

def targetP(accesionwLocation, accesionFinal, modification, MTS_targetP):

    # check wheter MTS or not!
    entrylist_targetP = list(MTS_targetP['Accession'])
    entryNumberList_targetP = list(MTS_targetP['Location'])

    # P49189 [275-293]
    # Q00610 [912-923]; P53675 [912-923]

    j = 0
    for entry in entrylist_targetP:
        if accesionFinal == entry:
            break
        else:
            j += 1

    entryNumber = int(entryNumberList_targetP[j-1])

    # locationLast = int(accesionwLocation.split('[')[1].split('-')[1].split(']')[0])
    locationFirst = int(accesionwLocation.split('[')[1].split('-')[0])

    if locationFirst <= entryNumber:
        try:
            url = f'https://www.ebi.ac.uk/proteins/api/proteins/{accesionFinal}'
            req = requests.get(url)
            result = req.json()
            geneSymbol = result['gene'][0]['name']['value']
        except:
            try:
                url = f'https://www.ebi.ac.uk/proteins/api/proteins/{accesionFinal}'
                req = requests.get(url)
                result = req.json()
                geneSymbol = result['gene'][0]['name']['value']
            except:
                geneSymbol = accesionFinal

        geneSymbolModification = f"{geneSymbol} [{accesionwLocation.split('[')[1]}"

        try:
            geneSymbolModification = geneSymbolModification.split(';')[0]
        except:
            pass

        targetPDecision = 'Yes'

        cache = {
            "Accession": accesionFinal,
            "Gene Symbol": geneSymbol,
            "Modifications": modification,
            "Positions in Master Proteins": geneSymbolModification,
            "Uniprot": 'No',
            "Uniprot Location": 0,
            "TargetP": targetPDecision,
            "TargetP Location": entryNumber,
        }

        return cache

def calculations(result, conditions, pairs):
    # ['WT', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO', 'KO']
    # KO/WT list2/list1
    testType = 'unpaired' # for only supports unpaired t test
    allAccessions = result['Accession']

    for pair in pairs:

        second = [i for i, x in enumerate(result.columns) if x == pair[0]] # KO list2
        first = [i for i, x in enumerate(result.columns) if x == pair[1]]  # WT list1

        num = 0

        for accesion in allAccessions:
            list1 = []
            for i in first:
                list1.append(result.iloc[:, i][num]) # to get all WT numbers in list

            list2 = []
            for i in second:
                list2.append(result.iloc[:, i][num]) # to get all KO numbers in list


            result.loc[num, f'Log2({pair[0]}/{pair[1]})'] = math.log(float(statistics.mean(list2)) / float(statistics.mean(list1)), 2)

            if testType == 'unpaired':
                result.loc[num, f'pvalue({pair[0]}/{pair[1]})'] = float(stats.ttest_ind(list2, list1, equal_var=True)[1])
                result.loc[num, f'-Log10 pvalue({pair[0]}/{pair[1]})'] = -math.log(float(stats.ttest_ind(list2, list1, equal_var=True)[1]),
                                                                      10)
            else:
                result.loc[num, f'-Log10 pvalue({pair[0]}/{pair[1]})'] = -math.log(float(stats.ttest_rel(list2, list1)[1]), 10)

            num += 1

    return result