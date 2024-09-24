
python3 ../SWEET-main/1.SWEET_sample_weight_calculating.py -f CLHLS_2012.Biomarkers.txt -s sample_weight.txt 
python3 ../SWEET-main/2.SWEET_edge_score_calculating.py -f CLHLS_2012.Biomarkers.txt -w sample_weight.txt -p individuals.txt -g Biomarkers.txt -s Individual/

#python3 ../SWEET-main/1.SWEET_sample_weight_calculating.py -f CLHLS_2012.Biomarkers.ref.txt -s sample_weight.txt
#python3 ../SWEET-main/2.SWEET_edge_score_calculating.py -f CLHLS_2012.Biomarkers.ref.txt -w sample_weight.txt -p Ref.individuals.txt -g Biomarkers.txt -s Individual/

#python3 ../SWEET-main/1.SWEET_sample_weight_calculating.py -f CLHLS_2012.Biomarkers.case.txt -s sample_weight.txt
#python3 ../SWEET-main/2.SWEET_edge_score_calculating.py -f CLHLS_2012.Biomarkers.case.txt -w sample_weight.txt -p Case.individuals.txt -g Biomarkers.txt -s Individual/
