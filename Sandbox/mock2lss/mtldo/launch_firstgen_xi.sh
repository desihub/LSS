for ty in ELG QSO LRG
do
	for i in {0..24}
	do
		python xifirstgen.py --type $ty --nran 10 --id $i
	done
done
