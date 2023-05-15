author: Manfred Scheucher, scheucher@math.tu-berlin.de

This repository provides sagemath scripts to visualize planar and plane graphs. 

=== draw.sage ===
Given a plain text file where each line encodes a planar graph, run "sage draw.sage [file] [informat] [outformat]" to compute a nice visualization for each graph. The script iteratively computes weighted Tutte embeddings and sets appropriate weights so that faces and edges become "nice". For more information checkout http://arxiv.org/abs/1708.06449. Possible input formats are "list" for edge list (python format), "s6" for sparse6 and "g6" for graph6. Possible output format are "pdf", "png", or "ipe". The latter is an XML format for the Ipe extensible drawing editor, which can be used to further modify the visualization in a "What You See Is What You Get" style. We refer to the official website https://ipe.otfried.org/

=== draw_plane.sage ===
When enumerating all embeddings of plane graphs using plantri, for example by the command
"plantri 5 -a -c1 -p > plane5.ascii", each line encodes a plane graph. More specifically, the cyclic rotations around each vertex are decoded. Given a plain text file where each line encodes a plane graph, run "sage draw_plane.sage [file]" to compute a visualization for each graph. The script adds auxiliary edges to triangulate the plane graph so that the embedding is unique (Whitney's theorem) and then computes coordinates of the vertices with tools for planar graphs. 
