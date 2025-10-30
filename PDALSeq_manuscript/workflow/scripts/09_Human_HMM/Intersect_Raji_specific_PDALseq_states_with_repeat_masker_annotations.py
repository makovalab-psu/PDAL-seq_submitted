import subprocess
import sys
import re

def run_enrichment_script(state_bed, annotation_bed, output_file, subsamples, simulations):
    # Construct the command to run the bash script
    command = f"bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {state_bed} {annotation_bed} {output_file} {subsamples} {simulations}"
    
    try:
        # Execute the command
        subprocess.run(command, shell=True, check=True)
        print(f"Successfully executed: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing the script: {e}")

def run_bed_generate_script(repeat_masker, repeat_class, output_file):
    # Construct the command to run the bash script
    command = f"bash workflow/scripts/07_PhyloHGMP_model/bed_generate.sh {repeat_masker} {repeat_class} > {output_file}"
    
    try:
        # Execute the command
        subprocess.run(command, shell=True, check=True)
        print(f"Successfully executed: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing the script: {e}")


Class = ["DNA","DNA?","DNA/Crypton","DNA/Crypton-A","DNA/hAT","DNA/hAT?","DNA/hAT-Ac","DNA/hAT-Blackjack","DNA/hAT-Charlie","DNA/hAT-hAT19","DNA/hAT-Tag1","DNA/hAT-Tip100","DNA/hAT-Tip100?","DNA?/hAT-Tip100?","DNA/Kolobok","DNA/Merlin","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/PiggyBac","DNA/TcMar","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc2","DNA/TcMar-Tigger","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1","LINE/L1-Tx1","LINE/L2","LINE/Penelope","LINE/RTE-BovB","LINE/RTE-X","Low_complexity","LTR","LTR/ERV1","LTR/ERVK","LTR/ERVL","LTR/ERVL-MaLR","LTR/Gypsy","RC/Helitron","Retroposon/SVA","rRNA","Satellite","Satellite/acro","Satellite/centr","Satellite/subtelo","scRNA","Simple_repeat","SINE/5S-Deu-L2","SINE/Alu","SINE/MIR","SINE/tRNA","SINE/tRNA-Deu","SINE/tRNA-RTE","snRNA","srpRNA","tRNA","Unknown"]

State_files = ["State_7_all", "State_7_high_PDAL-seq_expected", "State_7_variable_PDAL-seq_expected", "State_7_low_PDAL-seq_expected", "State_7_nonsyntenic", "State_7_syntenic"]

y = "State_7_low_PDAL-seq_expected"
x = "DNA"

for y in State_files:
	query_bed = "data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/"+y+".bed"
	rm_file = "resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out"
	subprocess.run("mkdir data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/repeat_masker_intersects/"+y, shell=True, check=True)
	out_bed = "data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/repeat_masker_beds/"+re.sub("/", "-", x)+".bed"
	for x in Class:
		CLASS = x
		out_enrichment = "data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/repeat_masker_intersects/"+y+"/"+re.sub("/", "-", x)+".csv"
		run_bed_generate_script(rm_file, CLASS, out_bed)
		run_enrichment_script(query_bed, out_bed, out_enrichment, 10000, 100)





