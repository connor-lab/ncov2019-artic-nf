#!/usr/bin/env python

import os, sys, time

def files_to_timestamp(path):
    files = [os.path.join(path, f) for f in os.listdir(path)]
    return dict ([(f, os.path.getmtime(f)) for f in files])

if __name__ == "__main__":

    path_to_watch = "/home/jyotirmoy/Documents/gms-artic/data/"
    #sys.argv[1]
    print('Watching {}..'.format(path_to_watch))

    before = files_to_timestamp(path_to_watch)

    while 1:
        time.sleep (2)
        after = files_to_timestamp(path_to_watch)

        added = [f for f in after.keys() if not f in before.keys()]
        removed = [f for f in before.keys() if not f in after.keys()]
        modified = []

        for f in before.keys():
            if not f in removed:
                if os.path.getmtime(f) != before.get(f):
                    modified.append(f)

        if added: os.system('./runner.sh')
        if removed: print('Removed: {}'.format(', '.join(removed)))
        if modified: os.system('./runner.sh')
            #print('Modified: {}'.format(', '.join(modified)))
            

        before = after