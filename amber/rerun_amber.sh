samp_name=$1

rm $(wc -l AMBER/${samp_name}*_ambermin.pdb | grep " 0 " | tr -s ' ' | cut -f 3 -d ' ')
cat fn_${samp_name}_premin.txt | xargs -n 1 -P 1 sh -c 'if [[ ! -f AMBER/${1%.pdb}_ambermin.pdb ]]; then /work/01872/nclement/scripts/amber/runAmber_single.sh $1 200; fi' sh
