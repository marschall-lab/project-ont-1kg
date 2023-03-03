#!/usr/bin/env python

import sys
import json

def bubble_chain_to_node_roles(bubble_chains):
    node_roles = defaultdict(set)
    for key in sorted(bubble_chains.keys()):
        chain = bubble_chains[key]
        if 'parent_chain' in chain:
            continue
        for node in chain['ends']:
            node_roles[node].add('chain_end')
        for bubble in chain['bubbles']:
            for node in bubble['ends']:
                node_roles[node].add('bubble_end')
            for node in bubble['inside']:
                node_roles[node].add('inside')
    return node_roles

def node_roles_to_colors(node_roles):
    d = dict()
    for node, node_role in node_roles.items():
        if 'inside' in node_role:
            d[node] = 'blue'
        if 'bubble_end' in node_role:
            d[node] = 'orange'
        if 'chain_end' in node_role:
            d[node] = 'red'
    return d

print('Not intended to be run, just to store some code fragments to read BubbleGun output', file=sys.stderr)
sys.exit(1)

bubble_chains = json.load(open('chm13-90c.r518.noseq.bubbles.json'))
node_roles = bubble_chain_to_node_roles(bubble_chains)
node_colors = node_roles_to_colors(node_roles)
