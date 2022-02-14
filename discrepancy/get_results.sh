# echo type N d best mean
# for file in halton.N*.d* prgn.N*.d* random.N*.d* naive.N*.d*; do
#   echo $file >&2
#   #first=$(grep -e "TA_improved_bardelta -iter" -e "Anzahl" $file | tr -d "\n" | sed $'s/+ /\\n/g')
#   #echo "$first"
#   grep -e "TA_improved_bardelta -iter" -e "Anzahl" $file | tr -d "\n" | sed $'s/+ /\\n/g' \
#     | sed 's/.*trials 10 pts.\(.*\).N\(.*\).d\(.*\).Walbest \(.*\) mean \(.*\) \+Anz.*/\1 \2 \3 \4 \5/'
#   #echo $first | sed 's/.*trials 10 pts.\(.*\).N\(.*\).d\(.*\).Walbest \(.*\) \+mean \(.*\)  Anz.*/\1 \2 \3 \4 \5/'
# done

file_string=$1

GREP="grep --text"

echo type N d best prog run finished
# naive 100000 50 4.896781e-02 delta  9 *
for file in $file_string; do
#for file in /scratch/01872/nclement/pts.2/out/*.d5.*; do
  N=$(echo $file | sed 's/.*.N\([0-9][^\.]*\).*/\1/');
  d=$(echo $file | sed 's/.*\.d\([0-9]\+\).*/\1/');
  run=$(echo $file | sed 's/.*\([0-9]\+\).out/\1/');
  type=$(echo $file | sed 's~.*d[0-9]\+\.\([^\.]\+\)\..*.out~\1~');
  if [[ "$type" == "prgNisan" ]]; then
    type=prgn.nisan
  elif [[ "$type" == "prgMT" ]]; then
    type=prgn
  fi
  N=100000
  #type=prgn.nisan
  >&2 echo "'$type' '$N' '$d' '$run'"

  if $GREP -q 'bardelta best' $file; then
    $GREP -e "bardelta.*Anzahl" $file |\
      sed -e "s/^/$type $N $d $run /" \
          -e 's/shrink bardelta best \([^ ]\+\)  mean.*/\1 bardelta /' \
          -e "s/$d \([0-9]\+\) \(.*\)/$d \2 \1 */";
  else
    num_finished=$($GREP -c 'bardelta.*Result' $file)
    #shrink bardelta Result 0.428052 at 94641
    max=$($GREP 'bardelta.*Result' $file | sed 's/.*Result \(.*\) at.*/\1/' | sort -nr | head -n 1)
    if [ -z "$max" ]; then
      max='NA'
    fi
    echo $type $N $d $max bardelta $run $num_finished
  fi

  if $GREP -q '^best.*Anzahl' $file; then
    $GREP -e "^best .*Anzahl" $file |\
      sed -e "s/^/$type $N $d $run /" \
          -e 's/best \([^ ]\+\)  mean.*/\1 delta /' \
          -e "s/$d \([0-9]\+\) \(.*\)/$d \2 \1 */";
  else
    #Result 0.428052 at 94641
    num_finished=$($GREP -c '^Result' $file)
    max=$($GREP '^Result' $file | sed 's/.*Result \(.*\) at.*/\1/' | sort -nr | head -n 1)
    if [ -z "$max" ]; then
      max='NA'
    fi
    echo $type $N $d $max delta $run $num_finished
  fi
done
