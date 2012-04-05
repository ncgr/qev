cat args.txt | xargs --max-procs 6 -n 7 -I myargs bash -c 'myargs'

