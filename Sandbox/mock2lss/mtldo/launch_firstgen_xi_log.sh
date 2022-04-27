for ty in ELG QSO LRG
do
	for i in {11..24}
	do
		for kind in lin log
		do
			python xifirstgen_v2.py --tracer $ty --mockrea $i --nran 20 --bin_type $kind
		done
	
	done
done
