from Bio import SeqIO,  Phylo, AlignIO
from numpy import where
import numpy as np
import pandas as pd
import itertools

# Use this cell to write your code
from copy import deepcopy
import re
import numpy as np
from math import exp, log
import math
# Use this cell for your code
from numpy import array, diag
import numpy as np
from scipy.linalg import expm, eig, inv
from math import exp, sqrt
import sys


class Vertex:
    def __init__(self,name):
        self.name = name
        self.parent = self
        self.children = []
        self.neighbors = []
        self.timesVisited = 0
        self.degree = 0
        self.inDegree = 0
        self.outDegree = 0
        self.newickLabel = ""
        self.sequence = ""        

class Tree:
    def __init__(self): 
        # here we collect the attributes a tree class requires
        self.vertices = []
        self.leaves = []
        self.newickLabel = ""
        # vertice map
        self.vmap = {}
        # edge lenghts
        self.edgeLengthMapForDirectedEdges = {}
        self.edgeLengthMapForUndirectedEdges = {}
        # vertice list
        self.verticesForPostOrderTraversal = []
        self.verticesForPreOrderTraversal = []
        self.newickLabelForLeavesSet_flag = False
        self.root = ""
        self.sequenceLength = 0
        self.pi_rho = [0]*4
        self.parsimonyScore = 0
        self.logLikelihoodScore = 0
        self.DNA_map = {"A":0,"C":1,"T":2,"G":3}
        self.transitionMatrices = {}
        self.conditionalLikelihoodVectors = {}
    def GetEdgeLength(self,parent, child):
      # this is self explanatory: we basically have the edge lenght as a distance mesure for instance
        #print(parent.name)
        #print(child.name)
        return (self.edgeLengthMapForDirectedEdges[(parent.name, child.name)])
    def ContainsVertex(self,v_name):
      # check if a given vertex is alreay in the vertex dict
        is_vertex_in_tree = v_name in self.vmap.keys()
        return (is_vertex_in_tree)
    def AddVertex(self,v_name):
      # add vertex to the vertex dict
      # in needs to be a class of Vertex
        v = Vertex(v_name)
        self.vmap[v_name] = v
        self.vertices.append(v)
    def GetVertex(self,v_name):
      return (self.vmap[v_name])
    def AddEdge(self,parent_name, child_name, edge_length):
      # adds an edge between parent and child with a given edge lenght
      # add parent to the vertex map
        p = self.vmap[parent_name]
      # add children to the vertex map
        c = self.vmap[child_name]
      # assign the parent to the child 
        c.parent = p
      # assing children to the parent 
        p.children.append(c)
      # increase the indegree of a child by one
        c.inDegree += 1
      # increase the outdegree of a parent by one
        p.outDegree += 1
      # increase the degree of a child by one (for undirected tree)
        c.degree += 1
      # increase the degree of a parent by one (for undirected tree)
        p.degree += 1
      # add the parent, child and the edge length to the edgeLengthMapForDirectedEdges dict       
        self.edgeLengthMapForDirectedEdges[(parent_name,child_name)] = edge_length      
    def ReadEdgeList(self,edgeListFileName):
      # this function basically reads the csv file and populates the data structure
        file = open(edgeListFileName,"r")
        for line in file:
            if len(line.split(",")) != 3:
                print (line)
            u_name, v_name, length = line.split(",")
            u_name = u_name.strip()
            v_name = v_name.strip()
            length = float(length.strip())
            # add the vertexes to the vertex map
            if not self.ContainsVertex(u_name):
                self.AddVertex(u_name)
            if not self.ContainsVertex(v_name):
                self.AddVertex(v_name)
            # adds edge between the vertexes
            self.AddEdge(u_name, v_name, length)
        file.close()
        self.SetLeaves()        
        self.SetRoot()
        self.SetVerticesForPostOrderTraversal()
    def PrintEdgeList(self):
        lst = []
        for edge in self.edgeLengthMapForDirectedEdges.keys():
            length = self.edgeLengthMapForDirectedEdges[edge]
            #print (edge[0],edge[1],length)  
            lst.append([edge[0],edge[1],length])
        return lst

    def SetLeaves(self):
      # the vertexes with the vertex degree 1 are the leves
        for v in self.vertices:
            if v.degree == 1:
                self.leaves.append(v)
    def SetRoot(self):
      # the vertex with the indegree 0 is root
        for v in self.vertices:
            if v.inDegree == 0:
                self.root = v              
    def SetNewickLabelForLeaves(self):
      # sets the vertex name as newick label
        if (self.newickLabelForLeavesSet_flag):
            pass
        else:
            self.newickLabelForLeavesSet_flag = True
            for v in self.leaves:
                v.newickLabel = v.name
    def ResetTimesVisited(self):
        for v in self.vertices:
            v.timesVisited = 0    
    def SimulateSeqEvol(model="JC"):
      if self.root.sequence == "":
        pass
    def ComputeNewickFormat(self):
        print ("Computing newick format")
        self.SetNewickLabelForLeaves()
        if len(self.verticesForPostOrderTraversal) == 0:
          # get the hidden vertices
            self.SetVerticesForPostOrderTraversal()        
        for v in self.verticesForPostOrderTraversal:
          # get both the children
            child_left, child_right = v.children
          # get the distance/edge lenght between the left child and the parent  
            length_left = self.GetEdgeLength(v,child_left)
          # get the distance/edge lenght between the right child and the parent
            length_right = self.GetEdgeLength(v,child_right)
          # write it in a newick format
            v.newickLabel = "(" + child_left.newickLabel + ":" + str(length_left)
            v.newickLabel += "," + child_right.newickLabel + ":" + str(length_right) + ")"
        self.newickLabel = v.newickLabel + ";"         
    def SetVerticesForPostOrderTraversal(self):
      # start from leaves and visit all the vertices such that children are visited before parents
        for v in self.leaves:
            self.verticesForPostOrderTraversal.append(v)
        # create a deep copy of the leaves list so that the original remains unchanged 
        # : any changes made to a copy of object do not reflect in the original object.
        verticesToVisit = deepcopy(self.leaves)        
        while len(verticesToVisit) > 0:
          # goes over the vertices and removes them iteratively from the "verticesToVisit" list
          # always remove a vertex from the end of the list
            c = verticesToVisit.pop()
          # p is the parent of the leave
            p = c.parent
          # increase the timesVisited variable by one           
            p.timesVisited += 1
          # if both the children of a parent are traversed, append the parent to "verticesForPostOrderTraversal" list
            if p.timesVisited == 2:
                self.verticesForPostOrderTraversal.append(p)
                # if both the children are traversed and if the vertex is a hidden node (not the root),
                # add the hidden node to the of verticesToVisit list
                if p.inDegree == 1:
                    verticesToVisit.append(p)
    def SetVerticesForPreOrderTraversal(self):
      if len(self.verticesForPostOrderTraversal) == 0:
        self.SetVerticesForPostOrderTraversal()
      self.verticesForPreOrderTraversal = deepcopy(self.verticesForPostOrderTraversal)
      self.verticesForPreOrderTraversal.reverse()
    def AddSequences(self,sequences_dic):
      for v_name in sequences_dic.keys():
        v = self.GetVertex(v_name)
        v.sequence = sequences_dic[v.name]
      self.sequenceLength = len(sequences_dic[v.name])
    def ReadNewickFile(self,newickFileName):
        file = open(newickFileName,"r")
        newick_string = file.readline()
        file.close()
        hidden_vertex_ind = 1        
        # use escape character for matching parenthesis 
        rx = r'\([^()]+\)'
        numTries = 0        
        while "(" in newick_string:            
            numTries += 1
            # search for the paranthesis
            m = re.search(rx,newick_string)
            # returns a tuple containing all the subgroups of the match "()"
            string_match = m.group()            
            # remove ( and )
            siblings_string = string_match[1:-1]
            c_left_name_and_length, c_right_name_and_length = siblings_string.split(",")
            c_left_name, c_left_length = c_left_name_and_length.split(":")
            c_right_name, c_right_length = c_right_name_and_length.split(":")
            if not self.ContainsVertex(c_left_name):
                self.AddVertex(c_left_name)
            if not self.ContainsVertex(c_right_name):
                self.AddVertex(c_right_name)
            hidden_vertex_name = "h" + str(hidden_vertex_ind)
            self.AddVertex(hidden_vertex_name)            
            self.AddEdge(hidden_vertex_name, c_left_name, float(c_left_length))
            self.AddEdge(hidden_vertex_name, c_right_name, float(c_right_length))
            newick_string = newick_string.replace(string_match,hidden_vertex_name)
            hidden_vertex_ind += 1                                
        self.SetLeaves()
        self.SetRoot()
    def ComputeParsimonyScoreUsingFitchHartiganAlgorithm(self):
        self.parsimonyScore = 0
        # Iterate over sites
        for site in range(self.sequenceLength):
          # Initialize character sets for leaves
          CharSet_dic = {v.name : set([v.sequence[site]]) for v in self.leaves}        
          for p in self.verticesForPostOrderTraversal:
            # Iterate over hidden vertices
            if p.outDegree > 0:
              left_child, right_child = p.children
              CharSet_left_child = CharSet_dic[left_child.name] 
              CharSet_right_child = CharSet_dic[right_child.name]
              # if interesection of char sets for children is zero then char set for parent is union of char set of children
              if (len(CharSet_left_child.intersection(CharSet_right_child))==0):
                CharSet_dic[p.name] = CharSet_left_child.union(CharSet_right_child)
                self.parsimonyScore += 1
              # else char set for parent is intersection of char set of children
              else:
                CharSet_dic[p.name] = CharSet_left_child.intersection(CharSet_right_child)
        print ("Parsimony score is", self.parsimonyScore)
        return (None)
    def GetJCTransitionMatrixForJCModel(self,t):
      # transition matrix for Jukes-Cantor model
        P = np.zeros((4,4),dtype=float)
        for row in range(4):
          for col in range(4):
            if row==col:
                P[row][col] = 0.25 * (1 + 3*(exp(-1 * (t/0.75))))
            else:
                P[row][col] = 0.25 * (1 - exp(-1 * (t/0.75)))

        return (P)
    def SetJCTransitionMatrices(self):
      # store the JC transition matrix for each edge
      # iterate over the vertices
        for c in self.vertices:
          if c!=self.root:


            p = c.parent
            # obtain t
            if c!=p:
              t = self.GetEdgeLength(p,c)
              # get the JC transition matrix for each edge/t         
              P = self.GetJCTransitionMatrixForJCModel(t)
              # store it in the transitionMatrices dict: key--> child
              self.transitionMatrices[c.name] = P[:]
              #print(self.transitionMatrices)

    def SetPiRho(self,pi_rho):

        self.pi_rho = pi_rho[:]
    def ComputeLogLikelihood(self):           

        for site in range(self.sequenceLength):  
     
          for v in self.verticesForPostOrderTraversal:

              if v.outDegree == 0:
                  self.conditionalLikelihoodVectors[v.name] = [0.0]*4
                  self.conditionalLikelihoodVectors[v.name][self.DNA_map[v.sequence[site]]] = 1.0
   
              else: 
                                
                  c_left, c_right = v.children                  
                  P_left = self.transitionMatrices[c_left.name]
                  P_right = self.transitionMatrices[c_right.name] 

                  self.conditionalLikelihoodVectors[v.name] = [0.0]*4

                  for char_parent in range(4):
                      sum_left_char_parent = 0.0
                      sum_right_char_parent = 0.0
                      for char_child in range(4):

                          if math.isnan(P_left[char_parent][char_child]) == True:
                              P_left[char_parent][char_child]= 1.000
                          
                          if math.isnan(P_right[char_parent][char_child]) == True:
                              P_right[char_parent][char_child] = 1.000

                          sum_left_char_parent += P_left[char_parent][char_child] * self.conditionalLikelihoodVectors[c_left.name][char_child]
                          sum_right_char_parent += P_right[char_parent][char_child] * self.conditionalLikelihoodVectors[c_right.name][char_child]

                          
                          if math.isnan(sum_left_char_parent) == True:
                              sum_left_char_parent= 1.000
                          
                          if math.isnan(sum_right_char_parent) == True:
                              sum_right_char_parent = 1.000

                          self.conditionalLikelihoodVectors[v.name][char_parent] = (sum_left_char_parent)*(sum_right_char_parent)  


                          if math.isnan(self.conditionalLikelihoodVectors[v.name][char_parent]) == True:
                              sum_left_char_parent= 1.000
                          
                          if math.isnan(self.conditionalLikelihoodVectors[v.name][char_parent]) == True:
                              sum_right_char_parent = 1.000

                  min_value = min(self.conditionalLikelihoodVectors[v.name])

                  self.conditionalLikelihoodVectors[v.name] /= min_value 

                  if math.isinf(log(min_value)) != True:
                      self.logLikelihoodScore += log(min_value)  

          sum_root = 0.0

          for char_root in range(4):
              sum_root += self.pi_rho[char_root] * self.conditionalLikelihoodVectors[self.root.name][char_root]
          self.logLikelihoodScore += log(sum_root)


def NewickFileToedge_list(newickFile):
  T_3 = Tree()
  T_3.ReadNewickFile(newickFile)
  list_1 = T_3.PrintEdgeList()

  with open('edge_list.txt', 'w'): pass
  for t in list_1:
      T_3 = Tree()
      edge_list_file_name_T_3 = 'edge_list.txt'
      with open(edge_list_file_name_T_3, 'a') as f:
        f.write( t[0] + ", " + t[1] + ", " +  str(t[2]) + '\n')

def processing_fast(fasta_data):
  seqs = {}
  i =0
  for record in SeqIO.parse(open(fasta_data, "r"), "fasta"):
      seqs[record.name] = str(record.seq)


      seqs[record.name]=seqs[record.name].replace("\n","")
      seqs[record.name]=seqs[record.name].replace("N","A")
      seqs[record.name]=seqs[record.name].replace("R","A")
      seqs[record.name]=seqs[record.name].replace("Y","A")
      seqs[record.name]=seqs[record.name].replace("K","A")
      seqs[record.name]=seqs[record.name].replace("S","A")
      seqs[record.name]=seqs[record.name].replace("W","A")
      seqs[record.name]=seqs[record.name].replace("M","A")
      seqs[record.name]=seqs[record.name].replace("H","A")


  return seqs

def aic(loglikeihood , p):

  aic = -2 * loglikeihood + 2 * p
  return aic

def bic(loglikeihood , p, n):

  bic = (-2 * loglikeihood) + math.log(n) * p
  return bic


if __name__ == '__main__':

  edge_list = "edge_list.txt"
  fasta_data = sys.argv[1]
  newick = sys.argv[2]

  T_3 = Tree()
  edge_list_file_name_T_3 = edge_list
  fasta_data = fasta_data

  NewickFileToedge_list(newick)
  seqs = {}

  seqs = processing_fast(fasta_data)
  seqs_length = len(seqs[list(seqs.keys())[0]])
  p =0 
  matrix = "JC"
  if matrix == "JC":
    p=1
  elif matrix == "K2P":
    p = 2

  T_3.ReadEdgeList(edge_list_file_name_T_3)
  T_3.AddSequences(seqs)
  # Set model parameters
  T_3.SetJCTransitionMatrices()
  pi_rho = [0.25]*4
  T_3.SetPiRho(pi_rho)
  T_3.ComputeLogLikelihood()
  print ("Loglikelihood score is ", T_3.logLikelihoodScore)

  print("AIC is ", aic(T_3.logLikelihoodScore , p))
  print("BIC is ", bic(T_3.logLikelihoodScore , p, seqs_length))

