source ~/Physics/Package/Qiskit/Env313/bin/activate
python3.13 draw_graph.py input_graph.txt -o graph_000000.pdf \
	   --init grid --iter 1000 \
	   --pos-out pos_000000.txt \
	   --edge-length 1.0 --repulsion 1.2 --gravity 0.01 \
	   --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000001.txt -o graph_000001.pdf --pos-in pos_000000.txt --pos-out pos_000001.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000002.txt -o graph_000002.pdf --pos-in pos_000001.txt --pos-out pos_000002.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000003.txt -o graph_000003.pdf --pos-in pos_000002.txt --pos-out pos_000003.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000004.txt -o graph_000004.pdf --pos-in pos_000003.txt --pos-out pos_000004.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000005.txt -o graph_000005.pdf --pos-in pos_000004.txt --pos-out pos_000005.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000006.txt -o graph_000006.pdf --pos-in pos_000005.txt --pos-out pos_000006.txt --figsize 5x8 --node-size 0.30 --font-size 11
python3.13 draw_graph.py output_graph_000007.txt -o graph_000007.pdf --pos-in pos_000006.txt --pos-out pos_000007.txt --figsize 5x8 --node-size 0.30 --font-size 11
