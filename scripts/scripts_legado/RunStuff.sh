#!/bin/bash
# Use current working directory
#$ -cwd
#
# Join stdout and stderr
#$ -j y
#
# Run job through bash shell
#$ -S /bin/bash
#
#You can edit the scriptsince this line
#
# Your job name
#$ -o ErrorMessages/Test_job_v2$TASK_ID
#$ -e ErrorMessages/Test_job_v2$TASK_ID
#
# Send an email after the job has finished
#$ -m e
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
#
module load r/4.1.3

echo "Inicio del trabajo: $(date)"
start_time=$(date +%s)

NumberOfSNPsToTest=1000
rm SNP/LL$SGE_TASK_ID.txt

### See number of SNPs per file

NumberOfFiles=$( ls OrderedCADDDatasets*.txt | wc -l )

Sum[0]=0
SNPNumber[0]=0
Numbers[0]=0
for (( i=1; i <= $NumberOfFiles; i++ ))
do
SNPNumber[$i]=$( wc -l OrderedCADDDatasets$i.txt | awk '{print $1}' )
SNPNumber[$i]=$(( ${SNPNumber[$i-1]} + ${SNPNumber[$i]} ))
Numbers[$i]=$( echo "scale=3; ${SNPNumber[$i]} / $NumberOfSNPsToTest" | bc )
Residue[$i]=$( echo ${Numbers[$i]} | awk -F. '{print $2}' )
Integer=$( echo ${Numbers[$i]} | awk -F. '{print $1}' )
if [ ${Residue[$i]} -eq "000" ]
then
Numbers[$i]=$Integer
Sum[$i]=$(( ${Sum[$i-1]} + ${Numbers[$i]} ))
Flag[$i]=0
else
Numbers[$i]=$(( $Integer + 1 ))
Sum[$i]=$(( ${Sum[$i-1]} + ${Numbers[$i]} ))
Flag[$i]=1
fi
done

Begin=0
End=0

echo "SGE_TASK_ID = $SGE_TASK_ID"

for (( i = 1; i <= $NumberOfFiles; i++ ))
do

End=$(( $End + ${Numbers[$i]} ))

if [ $SGE_TASK_ID -le $End ]
then
FileNumber=$i
RowsToTake=$NumberOfSNPsToTest
if [ $SGE_TASK_ID -eq $End ]
then
RowsToTake=$(( ${Residue[$i]} - 1 ))
fi
break
fi

done



if [ $FileNumber -eq "1" ]
then
StartNumber=$(( ( $SGE_TASK_ID ) * $NumberOfSNPsToTest + 2 ))
else
StartNumber=$(( ( $SGE_TASK_ID - ${Sum[$FileNumber]} ) * $NumberOfSNPsToTest ))
fi

StartNumber=$(( ( $SGE_TASK_ID  - 1 ) * $NumberOfSNPsToTest + 2 - ( ${Sum[$FileNumber-1]} * $NumberOfSNPsToTest )  ))
EndNumber=$(( $StartNumber + $RowsToTake - 1 ))

echo "FileNumber = $FileNumber HeadNumber = $StartNumber TailNumber = $EndNumber"

for (( i = $StartNumber; i <= $EndNumber; i++ ))
do

Chromosome=$( head -n$i OrderedCADDDatasets$FileNumber.txt | tail -n1 | awk '{print $2}' )
Position=$( head -n$i OrderedCADDDatasets$FileNumber.txt | tail -n1 | awk '{print $3}' )

rsnumber=$( grep $Position /mnt/Timina/dortega/hlopezh/data/data.bim | grep "^$Chromosome\s" | awk '{print $2}' )

grep $rsnumber /mnt/Timina/dortega/hlopezh/data/frq_v4.txt > SNP/SNPData$SGE_TASK_ID.txt

### Script to calculate LL

Rscript --vanilla LikelihoodCalculations_August25_2022.R $SGE_TASK_ID

# Rscript here



done



end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Fin del trabajo: $(date)"
echo "Duración del trabajo: $execution_time segundos"