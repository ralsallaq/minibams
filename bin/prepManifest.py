#!/usr/bin/env python
import numpy as np
import pandas as pd
import argparse
import sys
import os


def removeWhiteSpace(row):
    for i, v in row.iteritems():
        row.loc[i] = v.strip() if isinstance(v, str) else v
    return row
    

def main():
    """fixes taxids to primary taxids and saves a seqInfo file"""
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--events_capture_file','-i', help='text file with the event and number of records found for abnormal event', required=True)
    args_parser.add_argument('--list_of_bams_file','-b', help='text file with the unique list of bams found for the sample from records', required=True)
    args_parser.add_argument('--events_locii_fixed','-o', help='text file with fixed multiplicity', required=True)
    args = args_parser.parse_args()

    
    events_capture_file = str(args.events_capture_file).strip()
    list_of_bams_file = args.list_of_bams_file
    events_locii_fixedF = args.events_locii_fixed

    print("inputs are ",events_capture_file, list_of_bams_file, events_locii_fixedF, file=sys.stderr)
    

# Will try to test for file existence via bash
#    found = False
#    for filename in os.listdir():
#        print(filename)
#        if not filename in [events_capture_file]:
#            continue
#        else:
#            print('file {} found'.format(filename))
#            found = True
#    if not found:
#        print('file {} not found!!'.format(events_capture_file))
#        sys.exit(1)

 
    try: 
        events_capture = pd.read_csv(events_capture_file, sep="\t", header=None)
        print("done reading {}\n".format(events_capture_file), file=sys.stderr) 
        print(events_capture, file=sys.stderr) 
    except:
        events_capture = pd.DataFrame()
        print("failed reading {}".format(events_capture_file), file=sys.stderr)

    print(events_capture, file=sys.stderr)

    assert(events_capture.columns.shape[0] == 7), "expected 7 columns in {} got {}".format(events_capture_file, events_capture.columns.shape[0])
    assert(events_capture.shape[0] == 1), "expected 1 entry in {} got {}".format(events_capture_file, events_capture.shape[0])
    events_capture.columns = ['dataDir','sampleBase','locus','type1','type2','eventFile','nRecords']
    print(events_capture, file=sys.stderr) #7 columns


    if events_capture.loc[0,'nRecords'] == 0:

        print("no records for the event found ... creating empty file", file=sys.stderr)
        fh = open(events_locii_fixedF, 'wt')
        fh.close()




    else:

        bams_list = pd.read_csv(list_of_bams_file, sep="\t", header=None)
        bams_list = bams_list.drop_duplicates()

        print("done reading {}".format(list_of_bams_file), file=sys.stderr) 
    
    
        assert(bams_list.columns.shape[0] == 3), "expected 3 columns in {} got {}".format(list_of_bams_file, bams_list.columns.shape[0])
        bams_list.columns = ['target','sample','bamDir']
        print(bams_list, file=sys.stderr) #3 columns
    
        bams_list.loc[:,:] = bams_list.apply(lambda r: removeWhiteSpace(r))
        events_capture.loc[:,:] = events_capture.apply(lambda r: removeWhiteSpace(r))
       
        print("there records for this event", file=sys.stderr)
        print(events_capture[['eventFile','nRecords']], file=sys.stderr)

        records = pd.read_csv("sample_records/"+events_capture.loc[0,'eventFile'], sep="\t", skip_blank_lines=True, names=range(100))
        idx_keep_cols = ~(records.isnull().sum(axis=0) == records.shape[0])
        # Skip rows with ChrA or OrientA
        idx_keep_rows = ((records == 'ChrA') | (records == 'OrientA')).sum(axis=1) == 0
        records = records.loc[idx_keep_rows, idx_keep_cols]
        records.columns = range(records.shape[1])
        print(records.head(), file=sys.stderr)
        print(records.shape, file=sys.stderr)
 
        #### keep only the first 10 columns
        records = records.iloc[:,:10] 

        ### edit to remove file names within the first column
        records.loc[:,0] = records[0].apply(lambda r: r.split("txt:")[-1])
        records.loc[:,:] = records.apply(lambda r: removeWhiteSpace(r))
        records = records.drop_duplicates()
        records.index = range(records.shape[0])
        print(records.head(), file=sys.stderr)

        dupl_events = pd.concat([events_capture]*records.shape[0], axis=0, ignore_index=True)

        assert(dupl_events.shape[0] == records.shape[0]), "failed to duplicates events capture on records"

        merged1 = pd.concat([dupl_events, records], axis=1, ignore_index=True)

        print(merged1, file=sys.stderr)

        dupl_bams = pd.concat([bams_list]*merged1.shape[0], axis=0, ignore_index=True)
        dupl_merged1 = pd.concat([merged1]*bams_list.shape[0], axis=0, ignore_index=True)

        print(dupl_bams, file=sys.stderr)
        print(dupl_merged1, file=sys.stderr)

    
        events_locii = pd.concat([dupl_merged1, dupl_bams], axis=1, ignore_index=True)
        print(events_locii, file=sys.stderr)

        #assert(events_locii.drop_duplicates().shape[0] == events_locii.shape[0]), "unaccounted for duplication is present in events"

        assert(events_locii.columns.shape[0] == 20),'unexpected number of columns in manifest; expected 20 got {} columns'.format(events_locii.columns.shape[0])

        no_bams_idx = events_locii[19].isnull() 

        print('there are {} entries with no bams those might have the bams in column 15'.format(no_bams_idx.sum()), file=sys.stderr)

    
        if no_bams_idx.sum() > 0:
            assert(events_locii.loc[no_bams_idx,15].str.contains('.bam').sum() == no_bams_idx.sum()),'not all missing bams have bams in column 15, check if the bams are in another column'
            events_locii.loc[no_bams_idx,19] = events_locii.loc[no_bams_idx,15]
        
        print('after this step any null bams would make problems', file=sys.stderr)
        
        
        print(events_locii[17].isin(['WHOLE_GENOME', 'TRANSCRIPTOME', 'EXOME']).sum(), file=sys.stderr)
        assert(events_locii[17].isin(['WHOLE_GENOME', 'TRANSCRIPTOME', 'EXOME']).sum() == events_locii.shape[0]),'unknown analysis types or manifest is not in expected format'
        
        #### the assumption here is that _D* _A* _M* and _E* are tumor type and the rest are germline
        events_locii.loc[:,'sampleType'] = np.where((events_locii[18].str.contains('_D') | events_locii[18].str.contains('_A') | events_locii[18].str.contains('_E') | events_locii[18].str.contains('_M')),'tumor','germline')
        print(events_locii.head(), file=sys.stderr)
        print(events_locii.tail(), file=sys.stderr)
        
        print('checking if the events were mapped to 5 each corresponding to (2WGS, 2WES and 1RNA) per event', file=sys.stderr)
        
        if events_locii.shape[0] % 5 > 0: #the number of events are not multiples of 5 
        
            # Now we check that each event appears 5 times (2WGS, 2WES and 1RNA)
            # grouping by the columns defining the event
            # trying to deal with duplicate events is ugly, I will include all and then deduplicate by bed regions and bam file
            events = events_locii.groupby([1,2,3,4,6,7,8,9,10,11,12,13,14,15,16], as_index=False)
            
            okay_events_idx = (events.count() ==5)['sampleType']
            
            okay_events = events.count().loc[okay_events_idx]
            
            notOkay_events_idx = (events.count() != 5)['sampleType']
            
            print('there are {} okay events appearing 5 times corresponding to the 5 possible analysis types/sample types (2WGS, 2WES, 1RNA) and there are {} not okay events'.format(okay_events_idx.sum(), notOkay_events_idx.sum()), file=sys.stderr)
            
            ### those not okay events could be multiple germline mappings (e.g. SJMB009_H and SJMB009_G)
            ### to get the not okay events
            notokay_events = events.count().loc[notOkay_events_idx] 
            
            print('do the not okay events come from particular samples?', file=sys.stderr)
            print('samples with no okay events\\n',notokay_events[1].drop_duplicates(), file=sys.stderr)
            not_okay_full = pd.DataFrame()
            for s in notokay_events[1].values:
                notokay_events_s = notokay_events.loc[notokay_events[1]==s,[1,2,3,4]].drop_duplicates().values[0]
                not_okay_full = not_okay_full.append(events_locii[(events_locii[1]==notokay_events_s[0]) & (events_locii[2]==notokay_events_s[1]) & (events_locii[3]==notokay_events_s[2]) &(events_locii[4]==notokay_events_s[3])])
            
            not_okay_full = not_okay_full.drop_duplicates()
            
            print('the following events do not appear exactly  5 times', file=sys.stderr)
            print(not_okay_full, file=sys.stderr)
            
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
            
            
            #### remove the not_okay_full from events_locii
            fulldata = events_locii.shape[0]
            toreomve = not_okay_full.shape[0]
            events_locii = pd.concat([events_locii,not_okay_full], axis=0).drop_duplicates(keep=False)
            
            assert(fulldata-toreomve == events_locii.shape[0]),'removal of not-okay data was not successful'
            
            events_locii = events_locii.append(deduplicated_events)
            
            print('shape of deduplicated data = ',events_locii.shape, file=sys.stderr)
        
        print('at this point we have a file with correct multiplicity of events', file=sys.stderr)
        
        events_locii.to_csv(events_locii_fixedF, sep='\t', index=False, header=None)

if __name__ == "__main__":
    main()
