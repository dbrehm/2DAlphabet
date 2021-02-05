# Requires that you've run a b-only fit via run_MLfit.py at least once! Need to use the output.

from TwoDAlphabetClass import TwoDAlphabet
import sys, traceback, os, ROOT
from optparse import OptionParser
import subprocess
import header

def runLimit(twoDs,postfitWorkspaceDir,blindData=True,location=''):
    # Set verbosity - chosen from first of configs
    verbose = ''
    if twoDs[0].verbosity != False:
        verbose = ' -v '+twoDs[0].verbosity

    # Set systematics
    syst_option = ''
    for twoD in twoDs:
        for proc in twoD.inputConfig['PROCESS'].keys():
            if type(twoD.inputConfig['PROCESS'][proc]) == dict:
                if twoD.inputConfig['PROCESS'][proc]['CODE'] == 0:
                    if len(twoD.inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0: # If at any point there's a process
                        syst_option = ''    

    # Set signal strength range
    sig_option = ' --rMin 0 --rMax 5'

    # Run blind (turns off data everywhere) but don't mask (examines signal region)
    if blindData:
        blind_option = ' --run blind'
    else:
        blind_option = ''

    # Set the project directory
    if len(twoDs) > 1:
        identifier = twoDs[0].tag
        projDir = twoDs[0].tag
    else:
        identifier = twoDs[0].name
        projDir = twoDs[0].projPath

    card_name = 'card_'+identifier+'.txt'
    # Check if we can import post-fit result made during MLfit step
    if not os.path.isfile(postfitWorkspaceDir+'/fitDiagnostics.root'):
        print 'ERROR: '+postfitWorkspaceDir+'/fitDiagnostics.root does not exist. Please check that run_MLfit.py finished correctly. Quitting...'
        quit()

    for t in twoDs:
        del t

    # Make a prefit workspace from the data card
    print 'cd '+projDir
    with header.cd(projDir):
        t2w_cmd = 'text2workspace.py -b '+card_name+' -o workspace.root --X-no-jmax' 
        header.executeCmd(t2w_cmd)
        # header.setSnapshot(os.environ['CMSSW_BASE']+'/src/2DAlphabet/'+postfitWorkspaceDir+'/')

    print 'Making new workspace with snapshot'
    # Morph workspace according to imported fit result
    prefit_file = ROOT.TFile(projDir+'/workspace.root','UPDATE')
    print 'Getting workspace w'
    postfit_w = prefit_file.Get('w')
    print 'Opening fitDiagnostics'
    fit_result_file = ROOT.TFile.Open(postfitWorkspaceDir+'/fitDiagnostics.root')
    print 'Getting b only fit result'
    fit_result = fit_result_file.Get("fit_b")
    print 'Making RooArgSet out of fit result parameters'
    postfit_vars = ROOT.RooArgSet(fit_result.floatParsFinal())
    print 'Saving snapshot'
    postfit_w.saveSnapshot('initialFit',postfit_vars,True)
    # print 'Writing '+projDir+'limitworkspace.root'
    # fout = TFile(projDir+'/limitworkspace.root',"recreate")
    print 'Writing'
    prefit_file.WriteTObject(postfit_w,'w','Overwrite')
    print 'Closing'
    prefit_file.Close()

    #for idx in range(postfit_vars.getSize()):
    #    par_name = postfit_vars[idx].GetName()
    #    if postfit_w.var(par_name):
    #        print 'Setting '+par_name+' to '+str(postfit_vars[idx].getValV())+' +/- '+str(postfit_vars[idx].getError())
    #        var = postfit_w.var(par_name)
    #        var.setVal(postfit_vars[idx].getValV())
    #        var.setError(postfit_vars[idx].getError())

    # prefit_file.Close()

    print 'Getting current dir'
    current_dir = os.getcwd()

    aL_cmd = 'combine -M AsymptoticLimits workspace.root --snapshotName initialFit --saveWorkspace' +blind_option + syst_option# + sig_option 

    # Run combine if not on condor
    if location == 'local':    
        print 'cd '+projDir
        with header.cd(projDir):
            header.executeCmd(aL_cmd)
    # If on condor, write a script to run (python will then release its memory usage)
    elif location == 'condor':
        # Do all of the project specific shell scripting here
        shell_finisher = open('shell_finisher.sh','w')
        shell_finisher.write('#!/bin/sh\n')
        shell_finisher.write('cd '+projDir+'\n')
        shell_finisher.write(aL_cmd+'\n')
        shell_finisher.write('cd '+current_dir+'\n')
        shell_finisher.write('tar -czvf '+identifier+'.tgz '+projDir+'/\n')
        shell_finisher.write('cp '+identifier+'.tgz $CMSSW_BASE/../')
        shell_finisher.close()


parser = OptionParser()

parser.add_option('-q', '--tag', metavar='F', type='string', action='store',
                default =   '',
                dest    =   'quicktag',
                help    =   'Assigns a tag for this run')
parser.add_option('-d', '--projDir', metavar='F', type='string', action='store',
                default =   '',
                dest    =   'projDir',
                help    =   'Points to the directory where the b-only fit result is located')
parser.add_option("--unblindData", action="store_true", 
                default =   False,
                dest    =   "unblindData",
                help    =   "Unblind the observation and calculate the observed limit")
parser.add_option("--recycleAll", action="store_true", 
                default =   False,
                dest    =   "recycleAll",
                help    =   "Recycle everything from the previous run with this tag")


(options, args) = parser.parse_args()

inputConfigsAndArgs = args

# Assign and summarize input
stringSwaps = {}
inputConfigs = []
postfitWorkspaceDir = options.projDir
for i,c in enumerate(inputConfigsAndArgs):
    if ':' in c: # specify string swap
        stringSwaps[c.split(':')[0]] = c.split(':')[1]
        print c.split(':')[0] +' = '+c.split(':')[1]
    elif '.json' in c:
        inputConfigs.append(c)


print 'tag                              = ' + str(options.quicktag)
print 'Location of b-only fit workspace = ' + postfitWorkspaceDir 
print 'Unblind data points                = '+ str(options.unblindData)
print 'recycleAll\t = '+str(options.recycleAll)
print 'Config Replacements:'
for s in stringSwaps.keys():
    print '\t'+ s + ' -> ' + stringSwaps[s]
print 'Configs:'
for c in inputConfigs:
    print '\t'+c
    
# Check for workspace to load
if os.path.isfile(postfitWorkspaceDir +'/fitDiagnostics.root'):
    print 'Loading '+postfitWorkspaceDir +'/fitDiagnostics.root'
else:
    print "ERROR: post-fit '"+postfitWorkspaceDir+"/fitDiagnostics.root' could not be found. Please specify the location of the post-fit workspace (postfitb_workspace.root) that you would like to load. Quitting..."
    quit()


# Initialize
twoDinstances = []

# If simultaneous fit
if len(inputConfigs) > 1:
    # Instantiate all class instances
    for i in inputConfigs:
        instance = TwoDAlphabet(i,options.quicktag,options.recycleAll,stringSwaps=stringSwaps)
        twoDinstances.append(instance)

    # For each instance, check tags match and if they don't, ask the user for one
    for t in twoDinstances:
        if t.tag != twoDinstances[0].tag:
            raise ValueErro('ERROR: tags in configuration files do not match. '+t.tag+' does not match '+twoDinstances[0].tag+'. Please make sure they match and try again. Quitting...')
            
    thistag = twoDinstances[0].tag

    # Combine the cards 
    print 'cd ' + thistag
    with header.cd(thistag):
        card_combination_command = 'combineCards.py --X-no-jmax'
        for i in twoDinstances:
            card_combination_command += ' '+i.name+'/card_'+i.name+'.txt'
        card_combination_command += ' > card_'+thistag+'.txt'

        print 'Executing ' + card_combination_command
        subprocess.call([card_combination_command],shell=True)
        for num in range(1,len(twoDinstances)+1):
            subprocess.call(["sed -i 's/ch"+str(num)+"_//g' card_"+thistag+".txt"],shell=True)

    if 'condor' in os.getcwd(): runLimit(twoDinstances,postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='condor')
    else: runLimit(twoDinstances,postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='local')

# If single fit
else:
    instance = TwoDAlphabet(inputConfigs[0],options.quicktag,options.recycleAll,stringSwaps=stringSwaps)
    if 'condor' in os.getcwd(): runLimit([instance],postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='condor')
    else: runLimit([instance],postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='local')
