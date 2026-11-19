echo "Project quota:"
lfs quota -hp `stat -c "%g" /cfs/klemming/projects/supr/snic2021-5-492` /cfs/klemming

echo "Personal quota:"
lfs quota -hp $UID /cfs/klemming
