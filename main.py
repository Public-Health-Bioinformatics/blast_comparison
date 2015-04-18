#!/bin/env python

##############################################
# CompareBLASTs                              #
# A tool to compare the found hits from two  #
# BLAST searches with the same search query. #
#                                            #
# by Philipp B. Rentzsch                     #
# BCCDC Vancouver, BC                        #
# October 2014 - January 2015                #
##############################################

from __future__ import print_function
from time import strptime # convert string into time object
import optparse # commandline parsing
from BlastHit import * # BlastHit.py file

def loadBlastHits(file):
    '''
    Read a tabular BLAST file into a list of BlastHits
    
    file = (string) filename of tabular blast result file
    '''

    blastfile = open(file).readlines()
    return [BlastHit(hit) for hit in blastfile]
    
class CompareBLASTs(object):
    
    def __init__(self, old_hits, new_hits, email, name):
        '''
        Initialize the comaprison object:
        
        old_hits = List of BlastHits from the older BLAST Search
        new_hits = List of BlastHits from the newer, second BLAST Search
        email = Your email address, needed for use of NCBIs Enterez to prevent misuse of their service
        name = Query name that lead to the BlastHits to identify them later
        '''
    
        self.input_old_hits = old_hits
        self.input_new_hits = new_hits
        
        self.email = email
        self.name = name
        
    def compare(self):
        '''
        compares the two lists of BlastHits for more or less similar elements and extracts those elements form both lists that have no companion in each other
        '''
        
        # compare for exact (or similar) hits
        self.new_hits, self.old_hits = compareBLASTs(self.input_new_hits, self.input_old_hits)
        
        # retrieve basic information of the coresponding genes for all old hits
        self.oldGeneIDs = getIdList(self.old_hits['all'], self.email)
        
         # check all the old hits without a copy in the new hit list what happend to their associated gene (whether it still exists, was updated (=replaced) or deleted (=suppressed)
        oldOnly = {'live':[], 'replaced':[], 'suppressed':[]} # a bit confusing: live and lost are handled here equivalent since a hit that is live (=still existing in the db) but not found in the new BLAST search was 'lost' at some point
        for hit in self.old_hits['unknown']:
            for ID in hit.ids:
                if ID.db == 'gi':
                    oldOnly[self.oldGeneIDs[ID.num]['Status']].append(hit)
                    hit.status = self.oldGeneIDs[ID.num]['Status']
                    break
        
        self.new_hits['replacement'] = [] # equivalent to old_hits 'replaced'
        self.old_hits['lost'] = oldOnly['live']
        self.old_hits['suppressed'] = oldOnly['suppressed']
        self.old_hits['replacement'] = []
        
        # check the old hits with a known replacement tag, wether a replacing hit can be found in the new hits
        for num, hit in enumerate(oldOnly['replaced']):
            for ID in hit.ids:
                if ID.db == 'gi':
                    new_id = self.oldGeneIDs[ID.num]['ReplacedBy']
            found = False
            for num2, hit2 in enumerate(self.new_hits['unknown']):
                if new_id in [ID.num for ID in hit2.ids]:
                    same, differences = hit.compareHit(hit2, check_ids=False)
                    if same:
                        rep = self.new_hits['unknown'].pop(num2)
                        rep.status = 'replacement'
                        self.new_hits['replacement'].append(rep)
                        self.old_hits['replacement'].append(oldOnly['replaced'][num])
                        found = True
                        break
            if not found:
                # hit can be replaced but the replacement was nevertheless not found in the new Blast Search => lost/live
                self.old_hits['lost'].append(oldOnly['replaced'][num])
                oldOnly['replaced'][num].status = 'live'
                
        # get the basic info for those hit in the new search, that have no know relative in the old search
        self.newGeneIDs = getIdList(self.new_hits['unknown'], self.email)
        
        # estimate the time of the old BLAST (or last used database update) search by looking for the creation of the youngest entree that match to the old hits
        date_oldsearch = max([strptime(record['CreateDate'], "%Y/%m/%d") for record in self.oldGeneIDs.values()])

        # check wether all new hits with no relative in the old Search are indeed new (there for created after the last of the old Hits). I never had this case but one can never know ...
        self.new_hits['new'] = []
        self.new_hits['old'] = []
        for hit in self.new_hits['unknown']:
            if strptime(self.newGeneIDs[hit.ids[0].num]['CreateDate'], "%Y/%m/%d") < date_oldsearch:
                self.new_hits['old'].append(hit)
                hit.status = 'strange'
            else:
                self.new_hits['new'].append(hit)
                hit.status = 'new'
    
    def outputComparison(self, output_types = [lambda x: print(x)], top = 0, long_output = False, adaptive = True):
        '''
        Prints (and or writes to a file) the output of the BLAST comparison
        
        output_types = List of output lambdas like 'lambda x: print(x)' and 'lambda x: output_file.write(''.join([x, '\n']))'
        top = The number of Hits (from the top score) that are of interest for the comparion (0 = all)
        long_output = longer, more describitive output
        adaptive = in adaptive mode only those categories are displayed that appear like if there are no new hits in the second BLAST, this is not dispalyed
        '''
        
        # determine the number of hits (in the interested interval) that belong to each category
        
        hits_per_category = {'equal' : 0, 'similar' : 0, 'live' : 0, 'replaced' : 0, 'suppressed' : 0, 'new' : 0, 'strange' : 0}
        
        if top == 0:
            top = min(len(self.old_hits['all']), len(self.new_hits['all']))
        
        for hit in self.old_hits['all'][:top]:
            hits_per_category[hit.status] += 1
        for hit in self.new_hits['all'][:top]:
            if hit.status in ['new', 'strange']:
                hits_per_category[hit.status] += 1
                     
                
        if long_output:
            category_names = {'equal' : 'Found in both BLASTs results:\t%i',
                        'similar' : 'Found in both BLASTs results with slight changes in eValue or description:\t%i',
                        'live' : 'Not showing up for unknown reasons in the second BLAST (probably low scores):\t%i',
                        'replaced' : 'Replaced/updated before the second BLAST:\t%i',
                        'suppressed' : 'Deleted/suppressed before the second BLAST:\t%i',
                        'new' : 'New hits added to the database for the second BLAST:\t%i',
                        'strange' : 'Hits that do only appear in the second BLAST that should have appeared in the first:\t%i'}
        else:
            category_names = {'equal' : 'Equal Hits:\t%i',
                        'similar' : 'Changed Hits\t%i',
                        'live' : 'Lost Hits\t%i',
                        'replaced' : 'Replaced Hits:\t%i',
                        'suppressed' : 'Deleted Hits:\t%i',
                        'new' : 'New Hits:\t%i',
                        'strange' : 'New appearing Hits:\t%i'}   
        
        # for the different output channels (write to file or print)
        for output in output_types:
            # Always print the query name as more than one query can be found in a single BLAST
            if self.name:
                output('Query:\t%s' % self.name)
                        
            if long_output:
                output('Total Hits from old search:\t%i' % len(self.old_hits['all']) )
                output('Total Hits from new search:\t%i' % len(self.new_hits['all']) )

                if top != len(self.old_hits['all']):
                    output('Among the top %i hits were:' %top)
                else:
                    output('From all hits were:')
            
            for key in ['equal', 'similar', 'live', 'replaced', 'suppressed', 'new', 'strange']:
                if not adaptive or hits_per_category[key] > 0: # in (default) adaptive mode, only those hit categories are displayed that appeare (example: if there is no replaced hit, the replaced hits column is not displayed
                    output(category_names[key] % hits_per_category[key])
                    
            # separate from following queries        
            output('\n')
            
    def export_hit_categories(self, categories, path = ''):
        '''
        Exports the given categories into files (format similar to the input .blast format with a status column added at the end)
        
        categories= string with comma ',' delimited categories (e.g: new,all_old to export all new Hits and all the hits from the old search)
        path = file path to the exported files
        '''
        
        categories = categories.split(',')
        for category in categories:
            hits = None
            if category == 'new':
                hits = self.new_hits['new']
            if category == 'equal':
                hits = self.old_hits['same']
            if category == 'similar':
                hits = self.old_hits['similar']
            if category == 'lost':
                hits = self.old_hits['lost']
            if category == 'replaced':
                hits = self.old_hits['replaced']
            if category == 'suppressed':
                hits = self.old_hits['suppressed']
            if category == 'all_old':
                hits = self.old_hits['all']
            if category == 'all_new':
                hits = self.new_hits['all']
            if hits:
                with open(path + self.name + '_' + category + '.blast', 'w+') as f: # the query name and category speciefies the file name (e.g. Query7_all_new.blast)
                    for hit in hits:
                        f.write(str(hit) + '\n')

def performComparison(opts):
    '''
    The main function that compares two BLAST files against the same Query-Sequence
    
    opts = parsed OptionsParser
    '''
    
    new_hits = {}
    old_hits = {}
    # load the hits from the two input files
    new_hits_all = loadBlastHits(opts.new_Blast)
    old_hits_all = loadBlastHits(opts.old_Blast)
    
    # sort all hits for their repective query (as one BLAST file can contain multiple queries
    for hit in new_hits_all:
        if hit.name in new_hits.keys():
            new_hits[hit.name].append(hit)
        else:
            new_hits[hit.name] = [hit]

    for hit in old_hits_all:
        if hit.name in old_hits.keys():
            old_hits[hit.name].append(hit)
        else:
            old_hits[hit.name] = [hit]
            
    # make sure that both files where against the same queries
    assert old_hits.keys() == new_hits.keys()
    
    # define how to output the (general) results (print to console and/or save to file)
    output_types = []
    if opts.verbose:
        output_types.append(lambda x: print(x))
    if opts.save_output:
        output_file = open(opts.output_path + opts.save_output, 'w+')
        output_types.append(lambda x: output_file.write(''.join([x, '\n']))) #somewhat complicated expression because file.write does not automatically add a line end character
    
    for key in old_hits.keys():
        blastComparison = CompareBLASTs(old_hits[key], new_hits[key], opts.email, key)
        blastComparison.compare()
        blastComparison.outputComparison(output_types, opts.top, opts.long_output, opts.adaptive)
        
        # export specified hit categories to file
        if opts.export:
            blastComparison.export_hit_categories(opts.export, path = opts.output_path)
        
    if opts.save_output:
        output_file.close()       

if __name__ == '__main__':
    
    # general description of the program
    usage = \
    '''
    %prog [options]
    Neccessary to provide are the two tabular BLAST files old (-o) and new (-n)
    '''

    op = optparse.OptionParser(usage = usage)
    op.add_option('-o', '--old', default=None, dest='old_Blast', help = 'the older tabular BLAST file (24 columns)')
    op.add_option('-n', '--new', default=None, dest = 'new_Blast', help = 'the newer BLAST file')
    
    op.add_option('-t', '--top', type='int', default=0, help='specify when only the top X (integer value) hits for each query are of interest')
    
    op.add_option('-v', '--verbose', action='store_true', dest='verbose', default=True, help='print everything')
    op.add_option('-q', '--quiet', action='store_false', dest='verbose', help='stay quiet')
    
    op.add_option('-s', '--save', default=None, dest='save_output', help='file where the output is saved')
    op.add_option('-p', '--put', default='', dest='output_path', help='the path where the saved output and/or exported hit files are stored')
    
    op.add_option('-l', '--longOutput', action='store_true', dest='long_output', default = False, help='enable long names in the output')
    
    op.add_option('-a', '--adaptive', action='store_true', dest='adaptive', default = True, help='only display those hit classes, that have elements')
    op.add_option('-A', '--notAdaptive', action='store_false', dest='adaptive', help='display all elements')
    
    op.add_option('-e', '--email', default='test@test.com', help='email address of the user to send him/her notice of excess use')
    
    op.add_option('-x', '--export', default=None, help='export specified hit categories (Example: "-x new,old_all,suppressed", Categories: "equal, similar, live, replaced, suppressed, new, strange, all_old and all_new)"')
    
    opts, args = op.parse_args()
    
    assert opts.old_Blast and opts.new_Blast
    
    # executes the analysing program
    performComparison(opts)
