#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: Automatically test branch-site models to input cDNA alignment and newick format tree
# Created by Xiangchen Li on 2018/10/27

import os
import sys
import time
import shutil
import traceback
import logging.handlers
from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count
from ete3 import EvolTree, TreeStyle, NodeStyle, faces, AttrFace


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog='branch_site_pipeline.py')
    parser.add_argument('-i', '--input', dest='input_file', action='store', default=None,
                        help='Input tree nodes file')
    parser.add_argument('-o', '--output', dest='output_dir', action='store', default='output',
                        help='Output directory [default output]')
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                        help="Give verbose output [default False]")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None,
                        help="Logfile location [default None]")
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
                        help="How many threads will be used? [default all]")
    return parser.parse_args()


def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    exc_info = ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    return exc_info


def load_parameters():
    descendant_dict = defaultdict()
    with open(input_file, 'r') as f:
        all_lines = f.readlines()
        aln_file = all_lines[0].strip()
        if not os.path.exists(aln_file):
            logger.error('Input cDNA alignment file does not exist: {0}'.format(aln_file))
            sys.exit(1)
        else:
            logger.info('Input cDNA alignment file: {0}'.format(aln_file))
        tree_file = all_lines[1].strip()
        if not os.path.exists(tree_file):
            logger.error('Input tree file does not exist', tree_file)
            sys.exit(1)
        else:
            logger.info('Input tree file: {0}'.format(tree_file))
        t = EvolTree(tree_file, format=0)
        for descendant in t.iter_descendants():
            descendant_dict[descendant.node_id] = str(descendant)
        root = t.get_tree_root()
        id_list = []
        for leaf in t.traverse('preorder'):
            id_list.append(leaf.node_id)
        select_nodes = []
        if len(all_lines) > 2:
            for each_line in all_lines[2:]:
                s = each_line.strip()
                if s:
                    select_nodes.append(s)
        if select_nodes:
            nodes_line = ', '.join(select_nodes)
            logger.info('Input nodes: {0}'.format(nodes_line))
            for node in select_nodes:
                if node not in t:
                    logger.error('Error node: {0}'.format(node))
                    sys.exit(1)
            if not t.check_monophyly(values=select_nodes, target_attr='name'):
                logger.error('Some nodes are not monophyletic, please check')
                sys.exit(1)
            common_ancestor = t.get_common_ancestor(select_nodes)
        else:
            common_ancestor = root
            logger.info('No specific node')
        run_list = []
        for s in common_ancestor.iter_descendants():
            run_list.append(s.node_id)
        logger.info('These node ids will be checked: {0}'.format(str(run_list)))
        return run_list, aln_file, tree_file, descendant_dict


def run_codeml(mark_id, aln_file, tree_file, sleep):
    logger.info('sub-process: {0}'.format(str(mark_id)))
    time.sleep(round(sleep / args.threads, 2))
    run_dir = os.path.join(output_dir, str(mark_id))
    os.makedirs(run_dir)
    tree = EvolTree(tree_file, format=0)
    tree.link_to_alignment(aln_file)
    tree.run_model('M0')
    tree.workdir = run_dir
    tree.mark_tree([mark_id], marks=['#1'])
    tree.run_model('bsA.' + str(mark_id))
    tree.run_model('bsA1.' + str(mark_id))
    ps = tree.get_most_likely('bsA.' + str(mark_id), 'bsA1.' + str(mark_id))
    rx = tree.get_most_likely('bsA1.' + str(mark_id), 'M0')
    bsA = tree.get_evol_model('bsA.' + str(mark_id))
    p_bsA = bsA.classes['proportions'][2]
    wfrg2a = bsA.classes['foreground w'][2]
    if ps < 0.05 and float(wfrg2a) > 1:
        result = [mark_id, ps, rx, p_bsA, 'positive selection']
    elif rx < 0.05 and ps >= 0.05:
        result = [mark_id, ps, rx, p_bsA, 'relaxation']
    else:
        result = [mark_id, ps, rx, p_bsA, 'no signal']
    return result


def layout(node):
    name_face = faces.AttrFace("name", fsize=10, fgcolor="#009000")
    # If node is a leaf, add the node name
    if node.is_leaf():
        faces.add_face_to_node(name_face, node, column=0)
    faces.add_face_to_node(AttrFace("node_id"), node, column=0)


def tree_layout(tree_file, ps_node_list):
    t = EvolTree(tree_file, format=0)
    style_other = NodeStyle()
    style_other['size'] = 6
    style_ps = NodeStyle()
    style_ps['fgcolor'] = '#ff0000'
    style_ps['size'] = 6
    for node in t.iter_descendants():
        descendant = t.get_descendant_by_node_id(node.node_id)
        if node.node_id in ps_node_list:
            descendant.img_style = style_ps
        else:
            descendant.img_style = style_other
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_leaf_name = False
    result_picture = os.path.join(output_dir, 'positive_selection_tree.png')
    t.render(result_picture, tree_style=ts)


def main():
    p = Pool(args.threads)
    (run_list, aln_file, tree_file, descendant_dict) = load_parameters()
    logger.info('Run all sub-processes in total {0} threads'.format(str(args.threads)))
    all_results = []
    for i in run_list:
        result = p.apply_async(run_codeml, args=[i, aln_file, tree_file, run_list.index(i)])
        all_results.append(result)
    p.close()
    p.join()
    result_file = os.path.join(output_dir, 'all_results.txt')
    ps_node_list = []
    with open(result_file, 'w') as f:
        header = 'id: selection_p-value relaxation_p-value ps_proportion signal\n'
        f.write(header)
        for m in all_results:
            m_list = m.get()
            mark_id = m_list[0]
            ps = m_list[1]
            rx = m_list[2]
            proportion = m_list[3]
            selection = m_list[4]
            if selection == 'positive selection':
                ps_node_list.append(mark_id)
            descendant = descendant_dict[mark_id]
            result_line = '{0}: {1} {2} {3} {4}{5}\n\n'.format(str(mark_id),
                                                               str(ps),
                                                               str(rx),
                                                               proportion,
                                                               selection,
                                                               descendant)
            f.write(result_line)
    tree_layout(tree_file, ps_node_list)


if __name__ == '__main__':
    # Run as script
    # Parse command-line
    args = parse_cmdline()
    # Set up logging
    logger = logging.getLogger('branch_site_pipeline.py: %s' % time.asctime())
    t0 = time.time()
    input_file = args.input_file
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    # Do we need a log file?
    if args.logfile is not None:
        try:
            log_stream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(log_stream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.error("Could not open %s for logging", args.logfile)
            sys.exit(1)
    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    output_dir = args.output_dir
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    main()
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
