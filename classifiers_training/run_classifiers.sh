#!/bin/bash

# This script is used to both send jobs to slurm for the training of all classifiers, and log them.
# A log of all successful trainings is kept at the bottom of the file.

# HOW TO USE: 	Define 'training_commands', an array of the commands, with their respective time limits
# Format: 		"command,time_limit"
# Python scripts:	abstracts_train.py, TRI_train.py, MoR_train.py

# Example (see LOG below for the ones that have been run):
training_commands=(
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --get_final_model True --num_epochs 5,3:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --get_final_model True --num_epochs 5,3:00:00"
)

# Write absolute path where ExTRI2 folder is installed
ABS_PATH=''

# Loop over each training command and submit a job
for cmd_with_time in "${training_commands[@]}"
do

    # Separate the command and the time
    IFS=',' read -r cmd time_limit <<< "$cmd_with_time"

    # Extract a job name from the command for logging
    job_name=$(echo $cmd | grep -oP '(?<=python ).*(?=_train\.py)')
    
    # Write a temporary sbatch script
    cat > temp_sbatch.sh << EOT
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --error=${job_name}_%j.err
#SBATCH --output=${job_name}_%j.out
#SBATCH --ntasks=1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=48
#SBATCH --time=$time_limit

module load gcc/10.2.0 python/3.9.1 rocm/5.1.1

cd $ABS_PATH/ExTRI2/classifiers_training
unset PYTHONPATH
source training_env/bin/activate

# Trick used to get pytorch to detect the GPUs
export PYTHONPATH=$ABS_PATH/ExTRI2/classifiers_training/classifiers_training/lib/python3.9/site-packages:$ABS_PATH/.local/lib/python3.9/site-packages

# Run the training command
echo "RUN SCRIPT: $job_name"
$cmd
EOT

    # Submit the job
    sbatch temp_sbatch.sh
    echo "command: $cmd"

    # Remove the temporary sbatch script if not needed
    rm temp_sbatch.sh
done


: <<'LOG'
MODEL EXAMPLES
TRI
   "python TRI_train.py --pretrained_model BioLinkBERT,05:00:00"
 
MoR
   "python MoR_train.py --pretrained_model BioLinkBERT,05:00:00"

_____________________________________
TRAININGS AFTER 3RD VALIDATION
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --get_final_model True --num_epochs 5,3:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --get_final_model True --num_epochs 5,3:00:00"
   
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,8:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,8:00:00"
   "python TRI_train.py --pretrained_model BioLinkBERT --num_epochs 5,8:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --num_epochs 5,8:00:00"
   
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 8,15:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --num_epochs 8,15:00:00"
   "python TRI_train.py --pretrained_model BioLinkBERT --num_epochs 8,15:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --num_epochs 8,15:00:00" 
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 10,17:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --num_epochs 10,17:00:00" 
   "python TRI_train.py --pretrained_model BioLinkBERT --num_epochs 10,17:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --num_epochs 10,17:00:00" 
   
   "python TRI_train.py --pretrained_model BiomedNLP --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model BioBERT --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model BiomedNLP --span True --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model BioBERT --span True --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   
   "python TRI_train.py --pretrained_model BioBERT --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   
   "python TRI_train.py --pretrained_model BiomedNLP --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model BioBERT --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model sciBERT --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model BiomedNLP --num_epochs 5 ,07:00:00"
   "python MoR_train.py --pretrained_model BioBERT --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model sciBERT --num_epochs 5,07:00:00"

   "python MoR_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model BioBERT --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model sciBERT --span True --num_epochs 5,07:00:00"
   "python MoR_train.py --pretrained_model BiomedNLP --span True --num_epochs 5,07:00:00"
   "python TRI_train.py --pretrained_model BiomedNLP-large --span True --num_epochs 5,10:00:00"
   "python TRI_train.py --pretrained_model BiomedNLP-large --num_epochs 5,10:00:00"
   "python MoR_train.py --pretrained_model BiomedNLP-large --span True --num_epochs 5,10:00:00"
   "python MoR_train.py --pretrained_model BiomedNLP-large --num_epochs 5,10:00:00"


_____________________________________
TRAININGS AFTER 1_2 VALIDATION
   "python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,04:00:00"
   "python TRI_train.py --pretrained_model BioLinkBERT --num_epochs 5,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --num_epochs 5,04:00:00"

   "python TRI_train.py --pretrained_model BioLinkBERT --fold_validation False --span True --num_epochs 7,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --fold_validation False --span True --num_epochs 7,04:00:00"
   "python TRI_train.py --pretrained_model BioLinkBERT --fold_validation False --num_epochs 7,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --fold_validation False --num_epochs 7,04:00:00"
 

_____________________________________
TRAININGS AFTER 1st VALIDATION
   "python TRI_train.py --pretrained_model BioLinkBERT,04:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT,04:00:00"

   "python TRI_train.py --pretrained_model BioLinkBERT --num_epochs 5,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --num_epochs 5,05:00:00"

   "python TRI_train.py --pretrained_model BioLinkBERT --span True,04:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True,04:00:00"

   "python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,05:00:00"


___________________________________
TRAININGS BEFORE 1st VALIDATION

TRI
   "python TRI_train.py --pretrained_model BiomedNLP,05:00:00"
   "python TRI_train.py --pretrained_model BioLinkBERT,05:00:00"
   "python TRI_train.py --pretrained_model BioBERT,05:00:00"
   "python TRI_train.py --pretrained_model sciBERT,05:00:00"

MoR
   "python MoR_train.py --pretrained_model BiomedNLP,05:00:00"
   "python MoR_train.py --pretrained_model BioLinkBERT,05:00:00"
   "python MoR_train.py --pretrained_model BioBERT,05:00:00"
   "python MoR_train.py --pretrained_model sciBERT,05:00:00"
 

__________________________________
TRAININGS BEFORE PYTORCH LIGHTNING

Log of all the successful trainings done

    #Abstracts    
    "python abstracts_train.py  --pretrained_model distilbert \
                                --train_loss_weights 0.5 0.5 \
                                ,06:00:00"
    "python abstracts_train.py  --pretrained_model distilbert \
                                --train_loss_weights 0.04 0.96 \
                                ,06:00:00"
    "python abstracts_train.py  --pretrained_model distilbert \
                                --train_loss_weights 0.02 0.98 \
                                ,06:00:00"
    "python abstracts_train.py  --pretrained_model BiomedNLP\
                                ,06:00:00"
    "python abstracts_train.py  --pretrained_model BiomedNLP \
                                --train_loss_weights 0.02 0.98 \
                                ,06:00:00"
    
    #TRI
    "python TRI_train.py    --pretrained_model distilbert\
                            ,06:00:00"
    "python TRI_train.py    --pretrained_model BiomedNLP \
                            ,06:00:00"
    "python TRI_train.py    --pretrained_model BiomedNLP \
                            --masked True \
                            ,06:00:00"
    "python TRI_train.py    --pretrained_model BiomedNLP \
                            --added_tokens [TF] [TG] \
                            ,06:00:00"
    "python TRI_train.py    --pretrained_model BiomedNLP \
                            --added_tokens [TF] [TG] --masked True \
                            ,06:00:00"
    "python TRI_train.py    --pretrained_model BiomedNLP \ 
                            --added_tokens [TF] [TG] [G] --masked True \
                            ,06:00:00"

    # MoR
    "python MoR_train.py    --pretrained_model distilbert \
                            ,03:00:00"
    "python MoR_train.py    --pretrained_model BiomedNLP \
                            ,03:00:00"

LOG
