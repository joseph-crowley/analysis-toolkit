{
    gROOT->ProcessLine(".L loopers_cpp.so");
    std::string FILEDIR = "/ceph/cms/store/group/tttt/Worker/crowley/output/Analysis_TTJetRadiation/230223_tt_bkg_Cutbased";


    // Category Data_2018
    TChain *chData_2018 = new TChain("SkimTree");
    std::string sample_strData_2018("Data_2018");
    chData_2018->Add((FILEDIR + "/2018A/Run2018A_34_of_63.root").data());
    process_chain(chData_2018, sample_strData_2018);


}
