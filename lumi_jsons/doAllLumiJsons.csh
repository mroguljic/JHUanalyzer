#python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataB_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_B.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataB2_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_B2.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataC_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_C.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataD_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_D.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataE_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_E.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataF_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_F.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataG_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_G.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataH_bstar16.root -r run -l luminosityBlock -o 2016/lumi_data16_H.json

python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataB_bstar17.root -r run -l luminosityBlock -o 2017/lumi_data17_B.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataC_bstar17.root -r run -l luminosityBlock -o 2017/lumi_data17_C.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataD_bstar17.root -r run -l luminosityBlock -o 2017/lumi_data17_D.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataE_bstar17.root -r run -l luminosityBlock -o 2017/lumi_data17_E.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataF_bstar17.root -r run -l luminosityBlock -o 2017/lumi_data17_F.json

python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataA_bstar18.root -r run -l luminosityBlock -o 2018/lumi_data18_A.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataB_bstar18.root -r run -l luminosityBlock -o 2018/lumi_data18_B.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataC1_bstar18.root -r run -l luminosityBlock -o 2018/lumi_data18_C1.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataC2_bstar18.root -r run -l luminosityBlock -o 2018/lumi_data18_C2.json
python printLumiJson.py root://cmseos.fnal.gov//store/user/lcorcodi/bstar_nano/rootfiles/dataD_bstar18.root -r run -l luminosityBlock -o 2018/lumi_data18_D.json


#compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_B.json 2016/lumi_data16_B_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_B2.json 2016/lumi_data16_B2_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_C.json 2016/lumi_data16_C_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_D.json 2016/lumi_data16_D_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_E.json 2016/lumi_data16_E_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_F.json 2016/lumi_data16_F_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_G.json 2016/lumi_data16_G_final.json
compareJSON.py --and golden_jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 2016/lumi_data16_H.json 2016/lumi_data16_H_final.json

compareJSON.py --and golden_jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt 2017/lumi_data17_B.json 2017/lumi_data17_B_final.json
compareJSON.py --and golden_jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt 2017/lumi_data17_C.json 2017/lumi_data17_C_final.json
compareJSON.py --and golden_jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt 2017/lumi_data17_D.json 2017/lumi_data17_D_final.json
compareJSON.py --and golden_jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt 2017/lumi_data17_E.json 2017/lumi_data17_E_final.json
compareJSON.py --and golden_jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt 2017/lumi_data17_F.json 2017/lumi_data17_F_final.json

compareJSON.py --and golden_jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 2018/lumi_data18_A.json 2018/lumi_data18_A_final.json
compareJSON.py --and golden_jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 2018/lumi_data18_B.json 2018/lumi_data18_B_final.json
compareJSON.py --and golden_jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 2018/lumi_data18_C1.json 2018/lumi_data18_C1_final.json
compareJSON.py --and golden_jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 2018/lumi_data18_C2.json 2018/lumi_data18_C2_final.json
compareJSON.py --and golden_jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 2018/lumi_data18_D.json 2018/lumi_data18_D_final.json
