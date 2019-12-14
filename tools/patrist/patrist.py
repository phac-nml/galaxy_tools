import os
import sys
import argparse
from Bio import Phylo


def walk_up (tips, curnode, pathlen, cutoff):
    """
    Recursive function for traversing up a tree.
    """
    pathlen += curnode.branch_length
    if pathlen < cutoff:
        if curnode.is_terminal():
            tips.append( (curnode.name, pathlen) )
        else:
            for c in curnode.clades:
                tips = walk_up(tips, c, pathlen, cutoff)
    return tips

    
def walk_trunk (curnode, cutoff, parents):
    """
    Find all tips in the tree that are within a threshold distance
    of a reference tip.
    """
    # first go down to parent and up other branch
    tips = []
    pathlen = curnode.branch_length # 0.0184788
    p = parents[curnode]
    
    for c in p.clades:
        if c == curnode: continue
        if c.is_terminal():
            if pathlen + c.branch_length < cutoff:
                tips.append( (c.name, pathlen+c.branch_length) )
        else:
            tips.extend(walk_up([], c, pathlen, cutoff))
    
    # next walk down trunk until path length exceeds cutoff or hit root
    while parents.has_key(p):
        curnode = p
        pathlen += p.branch_length # + 0.0104047
        p = parents[curnode]
        if pathlen >= cutoff: break
        for c in p.clades:
            if c == curnode:
                continue
            if c.is_terminal():
                if pathlen + c.branch_length < cutoff: # + 0.0503079
                    tips.append( (c.name, pathlen+c.branch_length) )
            else:
                tips.extend(walk_up([], c, pathlen, cutoff))
    return tips

    

def find_short_edges(tree, cutoff, keep_ties=True, minimize=False, returnlist=False):
    """
    Find the shortest edge from the earliest sequence of a patient to a 
    any sequence from any other patient.
    minimize = keep only edge from earliest seq to the closest other seq
    keep_ties = [to be used in conjunction with minimize]
                report all edges with the same minimum distance
    """

    # generate dictionary of child->parent associations
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents.update({child: clade})

    tips = tree.get_terminals()
    res = {}
    for tip1 in tips:
        # find the shortest distance in sequences that "cluster" with this one
        min_dist = 99999.
        tip2 = []
        for tipname, dist in walk_trunk (tip1, cutoff, parents):
            if minimize and dist < min_dist:
                min_dist = dist
                tip2 = [[tipname, dist]]
            else:
                tip2.append([tipname, dist])

        t1 = tip1.name
        for t2, dist in tip2:
            # sort tip names in lexico order
            key = (t1, t2) if t1 < t2 else (t2, t1)
            if key in res:
                continue
            res.update({key: dist})
            if minimize and keep_ties:
                # output only one edge
                break

    if returnlist:
        reslist = []
        for key, dist in res.iteritems():
            reslist.append((key[0], key[1], dist))
        return reslist

    return res


def main():
    parser = argparse.ArgumentParser(
        description='Generate clusters of tips from a tree that have a path length within '
                    'a maximum distance of each other.'
    )
    parser.add_argument('tree', help='<input> file containing Newick tree string.')
    parser.add_argument('cutoff', type=float, help='Maximum patristic distance.')
    parser.add_argument('outfile', default=None, help='<output> file to write results in CSV format.')
    parser.add_argument('--minimize', help='Report no more than one nearest neighbour per tip.', action='store_true')
    parser.add_argument('--keep_ties', help='If more than one tip has the same patristic distance, '
                        'report all as nearest neighbours.', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite existing output file.', action='store_true')
    args = parser.parse_args()

    assert args.cutoff > 0, 'Cutoff %f must be greater than 0.' % (args.cutoff, )

    if os.path.exists(args.outfile) and not args.overwrite:
        print 'Output file', args.outfile, 'already exists, use --overwrite.'
        sys.exit()

    outfile = open(args.outfile, 'w')
    outfile.write('tree,tip1,tip2,dist,is.tie\n')

    trees = Phylo.parse(args.tree, 'newick')
    for treenum, tree in enumerate(trees):
        results = find_short_edges(tree, args.cutoff)
        for key, dist in results.iteritems():
            outfile.write('%d,%s,%s,%f\n' % (treenum, key[0], key[1], dist))

    outfile.close()

if __name__ == "__main__":
    main()
