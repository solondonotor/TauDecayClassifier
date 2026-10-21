import ROOT

##################################################################################################
# Load custom C-functions
##################################################################################################
#ROOT.gInterpreter.Declare('#include "my_analysis.cpp"')

processList = {
    "p8_ee_Ztautau_ecm91": {}, # USE THIS LINE WHEN YOU WANT TO RUN OVER THE COMPLETE SAMPLE
    #"events_164192792": {}, # USE THIS LINE TO ONLY RUN ABOUT THIS SUBSET OF THE WHOLE SAMPLE
    #"events_198604879": {},
}

includePaths = ["Analyzer_RTDM.h"]

inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA" # Whole sample
#inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Ztautau_ecm91" # Subset of the sample
#inputDir    = "/afs/cern.ch/user/a/aumarten/work/fcc/samples"
outputDir   = "result_RTDM"

##################################################################################################
# Main class of the Analysis
##################################################################################################
class RDFanalysis():
    output_var = []
    
    def analysers(dframe):
        ##################################################################################################
        # Set up MC_history collections
        ##################################################################################################
        dframe0 = (
            dframe
            #.Filter("0 < rdfentry_ && rdfentry_<= 100000000") # USE EVEN SMALLER SUBSET OF THE SAMPLE (FAST CODE TESTING)
            .Filter("0 < rdfentry_ && rdfentry_<= 500000")
            
            .Alias("Parents","Particle#0.index") # parents
            .Alias("Daughters","Particle#1.index") # daughters
            
            .Define("MC_taus","MCParticle::sel_pdgID(15,true)(Particle)")
            
            .Define("tauMinusPionicMode",  "FCCAnalyses::MCDecayMode::tauMinusPionicModes(Particle, Daughters)")
         
        )
        
        RDFanalysis.output_var.extend(["tauMinusPionicMode"])
        return dframe0
    
    def output():
        '''
        Output variables which will be saved to output root file.
        '''

        branch_list=RDFanalysis.output_var

        return branch_list
