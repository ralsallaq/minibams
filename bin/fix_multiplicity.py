import pandas as pd
import numpy as np
import argparse

def main():
    """fixes taxids to primary taxids and saves a seqInfo file"""
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--events_locii_raw','-i', help='text file with raw events', required=True)
    args_parser.add_argument('--events_locii_fixed','-o', help='text file with fixed multiplicity', required=True)
    args = args_parser.parse_args()

    events_locii_rawF = args.events_locii_raw
    events_locii_fixedF = args.events_locii_fixed

    dfn = pd.read_csv(events_locii_rawF, sep='\t', header=None)
    
    ##### fix the raw manifest
    print('here the assumption is that the manifest has 18 columns')  
    assert(dfn.columns.shape[0] == 18),'unexpected number of columns in manifest; expecting 18 columns'
    no_bams_idx = dfn[17].isnull()
    
    print('there are {} entries with no bams those might have the bams in column 15'.format(no_bams_idx.sum()))
    
    if no_bams_idx.sum() >0 :
        assert(dfn.loc[no_bams_idx,15].str.contains('.bam').sum() == no_bams_idx.sum()),'not all missing bams have bams in column 15, check if the bams are in another column'
        dfn.loc[no_bams_idx,17] = dfn.loc[no_bams_idx,15]
    
    print('after this step any null bams would make problems')
    
    def removeWhiteSpace(row):
        for i, v in row.iteritems():
            row.loc[i] = v.strip()
        return row
    
    dfn.loc[:,[1,2,3,4,5]] = dfn[[1,2,3,4,5]].apply(lambda r: removeWhiteSpace(r))
    
    ####splitting the manifest columns into sets to format column of index 17 (column 18)  into three columns
    dfn1 = dfn.iloc[:,0:17]
    
    dfn2 = dfn[17].apply(lambda r: pd.Series(r.strip().split(' '))) #these are analysis type/sample name/bam path
    assert(dfn2.shape[0]==dfn.shape[0]),'splitting the columns was not successful'
    
    dfn2.columns = [17,18,19]
    dfn = pd.concat([dfn1,dfn2], axis=1)
    
    assert(dfn[17].isin(['WHOLE_GENOME', 'TRANSCRIPTOME', 'EXOME']).sum() == dfn.shape[0]),'unknown analysis types or manifest is not in expected format'
    
    #### the assumption here is that _D* _A* _M* and _E* are tumor type and the rest are germline
    dfn.loc[:,'sampleType'] = np.where((dfn[18].str.contains('_D') | dfn[18].str.contains('_A') | dfn[18].str.contains('_E') | dfn[18].str.contains('_M')),'tumor','germline')
    print(dfn.head())
    print(dfn.tail())
    
    print('checking if the events were mapped to 5 each corresponding to (2WGS, 2WES and 1RNA) per event')
    
    if dfn.shape[0] % 5 > 0: #the number of events are not multiples of 5 
    
        ###### Now we check that each event appears 5 times (2WGS, 2WES and 1RNA)
        ### grouping by the columns defining the event
        #### trying to deal with duplicate events is ugly, I will include all and then deduplicate by bed regions and bam file
        events = dfn.groupby([1,2,3,4,6,7,8,9,10,11,12,13,14,15,16], as_index=False)
        
        okay_events_idx = (events.count() ==5)['sampleType']
        
        okay_events = events.count().loc[okay_events_idx]
        
        notOkay_events_idx = (events.count() != 5)['sampleType']
        
        print('there are {} okay events appearing 5 times corresponding to the 5 possible analysis types/sample types (2WGS, 2WES, 1RNA) and there are {} not okay events'.format(okay_events_idx.sum(), notOkay_events_idx.sum()))
        
        ### those not okay events could be multiple germline mappings (e.g. SJMB009_H and SJMB009_G)
        ### to get the not okay events
        notokay_events = events.count().loc[notOkay_events_idx] 
        
        print('do the not okay events come from particular samples?')
        print('samples with no okay events\\n',notokay_events[1].drop_duplicates())
        not_okay_full = pd.DataFrame()
        for s in notokay_events[1].values:
            notokay_events_s = notokay_events.loc[notokay_events[1]==s,[1,2,3,4]].drop_duplicates().values[0]
            not_okay_full = not_okay_full.append(dfn[(dfn[1]==notokay_events_s[0]) & (dfn[2]==notokay_events_s[1]) & (dfn[3]==notokay_events_s[2]) &(dfn[4]==notokay_events_s[3])])
        
        not_okay_full = not_okay_full.drop_duplicates()
        
        print('the following events do not appear exactly  5 times')
        print(not_okay_full)
        
        deduplicated_events = pd.DataFrame()
        if not_okay_full[not_okay_full[17]=='TRANSCRIPTOME'].shape[0] == notOkay_events_idx.sum(): ### this indicates transcriptome is not duplicated but the rest has duplicate either germline or tumor runs for WGS or WES
            groupByTarget = not_okay_full.groupby(17)
            for t,dft in groupByTarget:
                if t=='TRANSCRIPTOME':
                    deduplicated_events = deduplicated_events.append(dft)
                    continue
                
                gl_temp = dft[dft['sampleType']=='germline']
                mask_gl = (~gl_temp[18].str.contains('_H')) & (gl_temp[18].str.contains('_G')) 
                gl_temp = gl_temp[mask_gl]
                tu_temp = dft[dft['sampleType']=='tumor']
                idx_duplicated = tu_temp[range(0,18)].duplicated()
                tu_temp = tu_temp.loc[~idx_duplicated]
                temp = tu_temp.append(gl_temp) 
                deduplicated_events = deduplicated_events.append(temp)
            assert((deduplicated_events.groupby([17,'sampleType']).count()[19]==3).sum()==5),'deduplication of not okay events is unsuccessful'
        
        elif not_okay_full[not_okay_full[17]=='EXOME'].shape[0] == 2*notOkay_events_idx.sum(): ### this indicates exome and transcriptome are not duplicated but WGS is
            groupByTarget = not_okay_full.groupby(17)
            for t,dft in groupByTarget:
                if t=='TRANSCRIPTOME' or t=='EXOME':
                    deduplicated_events = deduplicated_events.append(dft)
                    continue
                
                gl_temp = dft[dft['sampleType']=='germline']
                mask_gl = (~gl_temp[18].str.contains('_H')) & (gl_temp[18].str.contains('_G')) 
                gl_temp = gl_temp[mask_gl]
                tu_temp = dft[dft['sampleType']=='tumor']
                idx_duplicated = tu_temp[range(0,18)].duplicated()
                tu_temp = tu_temp.loc[~idx_duplicated]
                temp = tu_temp.append(gl_temp) 
                deduplicated_events = deduplicated_events.append(temp)
            assert((deduplicated_events.groupby([17,'sampleType']).count()[19]==3).sum()==5),'deduplication of not okay events is unsuccessful'
        
        elif not_okay_full[not_okay_full[17]=='WHOLE_GENOME'].shape[0] == 2*notOkay_events_idx.sum(): ### this indicates exome, transcriptome and whole_genome are all not duplicated
            deduplicated_events = not_okay_full.copy()
            assert((deduplicated_events.groupby([17,'sampleType']).count()[19]==3).sum()==5),'deduplication of not okay events is unsuccessful'
        
        
        #### remove the not_okay_full from dfn
        fulldata = dfn.shape[0]
        toreomve = not_okay_full.shape[0]
        dfn = pd.concat([dfn,not_okay_full], axis=0).drop_duplicates(keep=False)
        
        assert(fulldata-toreomve == dfn.shape[0]),'removal of not-okay data was not successful'
        
        dfn = dfn.append(deduplicated_events)
        
        print('shape of deduplicated data = ',dfn.shape)
    
    print('at this point we have a file with correct multiplicity of events')
    
    dfn.to_csv(events_locii_fixedF, sep='\t', index=False, header=None)

if __name__ == "__main__":
    main()
