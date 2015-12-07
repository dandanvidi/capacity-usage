import csv, re

class Node(object):

    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent
        self.children = []

    def __str__(self):
        return self.name + '\n\t' + '\n\t'.join(self.children)
        
    def RemoveChild(self, name):
        if name not in self.children:
            raise KeyError('Node %s does not have a child named %s' % (self.name, name))
        self.children.remove(name)
    
    def InsertChild(self, name, index=-1):
        if index == -1:
            # the index -1 means that the new node should be inserted
            # at the end of the list of children
            self.children.append(name)
        elif index < -1:
            # the +1 is meant to fix the behaviour of 'insert' when using
            # negative indices
            self.children.insert(index + 1, name)
        else:
            self.children.insert(index, name)

    def MoveChild(self, name, index=-1):
        self.RemoveChild(name)
        self.InsertChild(name, index)        

class Tree(object):
    ROOT_NAME = 'root'

    def __init__(self):
        self.node_dict = {Tree.ROOT_NAME: Node(Tree.ROOT_NAME)}

    @staticmethod
    def FromTMS(fp, depth=-1):
        """
            TMS files are built like a BFS of the tree, therefore we use a stack
            to reconstruct the tree structure.
            
            If the parameter depth is set to a positive value, all nodes in that
            level will be leaves (i.e. without children). Leaves can have 
            non-unique names.
        """
        tree = Tree()
        current_branch = [Tree.ROOT_NAME] # the current list of ancestors (BFS stack)
        for j, row in enumerate(csv.reader(fp, delimiter='\t')):
            level = min([i for i in xrange(len(row)) if row[i]])
            name = row[level]
            while len(current_branch) > level + 1:
                current_branch.pop()
            parent = tree.GetNode(current_branch[-1])

            if depth > 0 and level >= depth:
                raise ValueError('Level (%d) exceeds the prescribed depth of the tree (%d)' % (level, depth))

            if depth > 0 and level == depth - 1: # this should be a leaf
                parent.children.append(name)
                continue

            if name in tree.node_dict:
#                print j, name, tree.node_dict
                raise KeyError('There is more than one non-leaf node with the name: ' + name)
            
            tree.node_dict[name] = Node(name, parent.name)
            parent.children.append(name)
            current_branch.append(name)
        return tree
    
    def ToTMS(self, fp):
         csv_writer = csv.writer(fp, delimiter='\t')
         bfs_stack = [(Tree.ROOT_NAME, -1)]
         while bfs_stack != []:
             node_name, level = bfs_stack.pop()
             
             # strip the name from anything in curly braces
             stripped_name = re.sub('{[\{\}]*}', '', node_name)
             
             if level >= 0:
                 csv_writer.writerow([None]*level + [stripped_name])
             node = self.GetNode(node_name)
             if node is None: # this is a leaf
                 continue
             bfs_stack += [(child_name, level + 1) for child_name in reversed(node.children)]

    def GetNode(self, name):
        return self.node_dict.get(name, None)

    def MoveNode(self, name, new_parent=None, index=0):
        """
            More a node to a new location in the tree (together with the entire
            subtree). If the argument 'new_parent' is None, the node will only
            move to a new index inside its current list of siblings.
        """
        node = self.GetNode(name)
        if node is None:
            raise KeyError('The node %s does not exist' % name)
        if new_parent is None or new_parent == node.parent:
            # the parent is the same, only the order of siblings should change
            self.GetNode(node.parent).MoveChild(name, index) 
        else:
            # we need to move this node to a new parent (together with its subtree)
            self.GetNode(node.parent).RemoveChild(name) 
            parent = self.GetNode(new_parent)
            if parent is None:
                raise Exception('Cannot find parent node %s' % parent)
            parent.InsertChild(name, index)

    def RemoveNode(self, name):
        """
            Remove a node an all its subtree
        """
        node = self.GetNode(name)
        if node is None:
            raise KeyError('The node %s does not exist' % name)
        
        # delete all the node's children first (recursively)
        for child in node.children:
            if child in self.node_dict:
                self.RemoveNode(child)
            else:
                self.RemoveLeaf(name, child)

        # get the list of the children of the parent node (i.e. the siblings)
        parent = self.GetNode(node.parent)
        if parent is None:
            raise Exception('Cannot find parent node %s' % node.parent)
        parent.RemoveChild(name)

        # delete the node from the dictionary
        del self.node_dict[name]
        
    def MoveLeaf(self, name, parent_name, index=0):
        """
            Move a leaf among its siblings. To move a leaf to a new node,
            remove it from one parent and insert it in another.
        """
        parent = self.GetNode(parent_name)
        if parent is None:
            raise KeyError('The parent node %s does not exist' % name)
        parent.MoveChild(name, index) 
        
    def RemoveLeaf(self, name, parent_name):
        """
            Remove a leaf from its parent's list of children
        """
        parent = self.GetNode(parent_name)
        if parent is None:
            raise KeyError('The parent node %s does not exist' % name)
        parent.RemoveChild(name) 

    def RemoveAllLeaves(self, name):
        for node in self.node_dict.values():
            if name in node.children:
                node.RemoveChild(name)

    def InsertNode(self, name, parent_name, index=-1):
        """
            Insert a new node under a specific parent
        """
        if self.GetNode(name) is not None:
            raise KeyError('A node with this name already exists: ' + name)
        parent = self.GetNode(parent_name)
        if parent is None:
            raise KeyError('The parent node %s does not exist' % parent_name)
        self.node_dict[name] = Node(name, parent_name)
        parent.InsertChild(name, index)

    def InsertLeaf(self, name, parent_name, index=-1):
        """
            Insert a new leaf under a specific parent
        """
        parent = self.GetNode(parent_name)
        if parent is None:
            raise KeyError('The parent node %s does not exist' % name)
        parent.InsertChild(name, index)

    def RenameNode(self, name, new_name):
        """
            Rename a node. This function makes sure that the parent and 
            all the children of this node are updated.
        """
        node = self.GetNode(name)
        if node is None:
            raise KeyError('The node %s does not exist' % name)
        siblings = self.GetNode(node.parent).children
        siblings[siblings.index(name)] = new_name
        node.name = new_name
        del self.node_dict[name]
        self.node_dict[new_name] = node
        
        # need to rename the "parent" of all the children
        for child in map(self.GetNode, node.children):
            if child is not None:
                child.parent = new_name
        
    def MergeNodeInto(self, name_from, name_into):
        node_from = self.GetNode(name_from)
        node_into = self.GetNode(name_into)
        if node_from is None:
            raise KeyError('The node %s does not exist' % name_from)
        if node_into is None:
            raise KeyError('The node %s does not exist' % name_into)

        for child in node_from.children:
            if child in self.node_dict:
                self.MoveNode(child, name_into)
            else:
                node_from.RemoveChild(child)
                node_into.InsertChild(child)
        
    def __str__(self):
        return self.GetNode(Tree.ROOT_NAME).__str__()

if __name__ == "__main__":
    tree = Tree.FromTMS(open('../data/KO_gene_hierarchy_general.tms', 'r'), 4)
#    tree.ToTMS(open('tmp.tms', 'w'))
#    print tree
#    tree.MoveNode('Human Diseases', None, -2)
#    print tree
#    print tree.GetNode('Human Diseases')
#    tree.RenameNode('Human Diseases', 'Diseases')
#    tree.RemoveNode('Diseases')
#    print tree
#    tree.InsertNode('root', 'Human Diseases')
#    print tree
#    tree.InsertNode('Cholera', 'Human Diseases')
#    tree.InsertNode('Metabolism', 'Cancer')
#    tree.InsertLeaf('Metabolism', 'K00000')
#    tree.InsertLeaf('Human Diseases', 'K00000')
#    tree.MoveNode('Cancer', 'Human Diseases', 0)
#    print tree.GetNode('Human Diseases')
#    print tree.GetNode('Metabolism')