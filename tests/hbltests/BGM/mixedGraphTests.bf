/*
	Simulate data from mixed network and reconstruct network from these data
	from the following graph:
	
	  [A]   [B]   [C]   (D)
	         | \ / |
	         | / \ |
	        [E]   (F)
	           \ /  \
	           (G)   |
	            |   /
	             \ /
	             (H)
	
	where [X] is a discrete node and (Y) is a continuous (Gaussian) node.
	
*/


// settings
ncases = 100;
BGM_CONTINUOUS_MISSING_VALUE = -666;
pr_censor = 0.2;


path = HYPHY_BASE_DIRECTORY+".."+DIRECTORY_SEPARATOR
	+"tests"+DIRECTORY_SEPARATOR+"hbltests"+DIRECTORY_SEPARATOR
	+"BGM"+DIRECTORY_SEPARATOR;

ExecuteAFile (path+"bayesgraph.ibf");


// note: declaring NodeType is not necessary but is here for documentation purposes

// NParentCombs should be automatically populated 

mixedBN = {};
mixedBN["A"] = {"Parents":{}, "NodeType":0, "Levels":2, "NParentCombs":1};
mixedBN["B"] = {"Parents":{}, "NodeType":0, "Levels":2, "NParentCombs":1};
mixedBN["C"] = {"Parents":{}, "NodeType":0, "Levels":2, "NParentCombs":1};
mixedBN["D"] = {"Parents":{}, "NodeType":1, "NParentCombs":1};
mixedBN["E"] = {"Parents":{{"B"},{"C"}}, "NodeType":0, "Levels":2, "NParentCombs":4};
mixedBN["F"] = {"Parents":{{"B"},{"C"}}, "NodeType":1, "NParentCombs":4};
mixedBN["G"] = {"Parents":{{"E"},{"F"}}, "NodeType":1, "NParentCombs":2};
mixedBN["H"] = {"Parents":{{"F"},{"G"}}, "NodeType":1, "NParentCombs":1};


adjmat = list2adjmat (mixedBN);


// Initialize conditional PDFs with arbitrary hyperparameters.
// assume imaginary sample size of 1000


(mixedBN["A"])["CPDFs"] = {{"Random({{500, 500}}, {\"PDF\":\"Dirichlet\"})"}};
(mixedBN["B"])["CPDFs"] = {{"Random({{500, 500}}, {\"PDF\":\"Dirichlet\"})"}};
(mixedBN["C"])["CPDFs"] = {{"Random({{500, 500}}, {\"PDF\":\"Dirichlet\"})"}};

(mixedBN["D"])["CPDFs"] = {{"cg_params({{0.0}}, 1, 1, {{1}})",}};

(mixedBN["E"])["CPDFs"] = {{"Random({{200,50}}, {\"PDF\":\"Dirichlet\"})",
							"Random({{50,200}}, {\"PDF\":\"Dirichlet\"})",
							"Random({{50,200}}, {\"PDF\":\"Dirichlet\"})",
							"Random({{200,50}}, {\"PDF\":\"Dirichlet\"})"}};



// conditional Gaussian (CG) node with 2 discrete parents
(mixedBN["F"])["CPDFs"] = {{
	"cg_params({{1.0}}, 1, 1, {{1}})",
	"cg_params({{-1.0}}, 1, 1, {{1}})",
	"cg_params({{-1.0}}, 1, 1, {{1}})",
	"cg_params({{1.0}}, 1, 1, {{1}})"
	}};

// CG node with 1 discrete and 1 CG parent
(mixedBN["G"])["CPDFs"] = {{
	"cg_params({{-1.0, 1.0}}, 1, 1, {{1,0},{0,1}})",
	"cg_params({{1.0, -1.0}}, 1, 1, {{1,0},{0,1}})"
	}};

// CG node with 2 CG parents
(mixedBN["H"])["CPDFs"] = {{
	"cg_params({{0, 1.0, -1.0}}, 1, 1, {{1,0,0},{0,1,0},{0,0,1}})"
}};


// adds Parameters entry
instantiate_CPDFs(mixedBN);


// note - eigendecomposition routine becomes unstable for large matrices,
// which are induced for large data sets
sim = simulate_data(mixedBN, ncases);


// now try to recover the network
// use empirical summary statistics for hyperparameters
nodes = {};
nodes["0"] = add_discrete_node ("A", 2, 0, 2);
nodes["1"] = add_discrete_node ("B", 2, 0, 2);
nodes["2"] = add_discrete_node ("C", 2, 0, 2);
nodes["3"] = add_gaussian_node ("D", 2, 0, 0.0, 0.1, 2);
nodes["4"] = add_discrete_node ("E", 2, 0, 2);
nodes["5"] = add_gaussian_node ("F", 2, 0, 0.0, 0.1, 2);
nodes["6"] = add_gaussian_node ("G", 2, 0, 0.0, 0.1, 2);
nodes["7"] = add_gaussian_node ("H", 2, 0, 0.0, 0.1, 2);

num_nodes = Abs(nodes);
names = Rows(nodes);

BayesianGraphicalModel my_bgm = (nodes);

// have to introduce constraints for mixed network or else this 
// can cause serious problems for order MCMC
constraints = {num_nodes, num_nodes};
for (parent = 0; parent < num_nodes; parent+=1) {
	for (child = 0; child < num_nodes; child+=1) {
		if ( child != parent 
			  && ( (nodes[parent])["NodeType"] == 1 && (nodes[child])["NodeType"] == 0 ) ) {
			constraints[parent][child] = -1;
		}
	}
}
setConstraints ("my_bgm", constraints);


attach_data ("my_bgm", sim, 0, 0, 0);

cache0 = get_node_score_cache("my_bgm");


fprintf (stdout, "\nRUNNING GRAPH-MCMC\n");

result0 = graph_MCMC ("my_bgm", 500000, 500000, 100, 0);


temp = check_edgelist (result0, adjmat, 0.8);
fprintf (stdout, "\tTrue positives = ", temp[0], "\n");
fprintf (stdout, "\tFalse negatives = ", temp[1], "\n");
fprintf (stdout, "\tFalse positives = ", temp[2], "\n");
fprintf (stdout, "\tTrue negatives = ", temp[3], "\n");

sens = temp[0]/(temp[0]+temp[1]);
spec = temp[3]/(temp[2]+temp[3]);

fprintf (stdout, "\tSensitivity (TP/TP+FN) = ", sens, "\n");
fprintf (stdout, "\tSpecificity (TN/TN+FP) = ", spec, "\n");

fprintf (stdout, "Specificity > 75% and specificity > 90% for cutoff = 0.8 ... ");
if (sens > 0.75 && spec > 0.9) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


display_MCMC_chain (result0);
write_edgelist(path+"mixed.edges", result0, num_nodes, 1);
mcmc_graph_to_dotfile(path+"mixed.dot", 0.6, result0, nodes);



fprintf (stdout, "\nRUNNING ORDER-MCMC\n");

result1 = order_MCMC ("my_bgm", 100000, 100000, 1000);


temp = check_edgelist (result1, adjmat, 0.8); /* TP, FN, FP, TN */
fprintf (stdout, "\tTrue positives = ", temp[0], "\n");
fprintf (stdout, "\tFalse negatives = ", temp[1], "\n");
fprintf (stdout, "\tFalse positives = ", temp[2], "\n");
fprintf (stdout, "\tTrue negatives = ", temp[3], "\n");

sens = temp[0]/(temp[0]+temp[1]);
spec = temp[3]/(temp[2]+temp[3]);

fprintf (stdout, "\tSensitivity (TP/TP+FN) = ", sens, "\n");
fprintf (stdout, "\tSpecificity (TN/TN+FP) = ", spec, "\n");

fprintf (stdout, "Specificity > 75% and specificity > 90% for cutoff = 0.8 ... ");
if (sens > 0.75 && spec > 0.9) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


write_edgelist(path+"mixed_order.edges", result1, num_nodes, 1);
mcmc_graph_to_dotfile(path+"mixed_order.dot", 0.6, result1, nodes);



csim = {Rows(sim), Columns(sim)};
for (i = 0; i < Rows(csim); i+=1) {
	for (j = 0; j < Columns(csim); j+=1) {
		if (Random(0,1) < pr_censor) {
			if ((nodes[j])["NodeType"] == 0) {
				csim[i][j] = -1;
			} else {		
				csim[i][j] = BGM_CONTINUOUS_MISSING_VALUE;
			}
		} else {
			csim[i][j] = sim[i][j];
		}
	}
}

//fprintf (stdout, csim, "\n");
fprintf (stdout, "\nTEST IMPUTATION WITH 20% MISSING DATA\n");


attach_data ("my_bgm", csim, 100, 1, 100);
cache1 = get_node_score_cache ("my_bgm");

result2 = graph_MCMC ("my_bgm", 500000, 500000, 100, 0);


temp = check_edgelist (result2, adjmat, 0.8);
fprintf (stdout, "\tTrue positives = ", temp[0], "\n");
fprintf (stdout, "\tFalse negatives = ", temp[1], "\n");
fprintf (stdout, "\tFalse positives = ", temp[2], "\n");
fprintf (stdout, "\tTrue negatives = ", temp[3], "\n");

sens = temp[0]/(temp[0]+temp[1]);
spec = temp[3]/(temp[2]+temp[3]);

fprintf (stdout, "\tSensitivity (TP/TP+FN) = ", sens, "\n");
fprintf (stdout, "\tSpecificity (TN/TN+FP) = ", spec, "\n");

fprintf (stdout, "Specificity > 70% and specificity > 90% for cutoff = 0.8 ... ");
if (sens > 0.7 && spec > 0.9) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


display_MCMC_chain (result2);
write_edgelist(path+"mixed_missing.edges", result2, num_nodes, 1);
mcmc_graph_to_dotfile(path+"mixed_missing.dot", 0.6, result2, nodes);

