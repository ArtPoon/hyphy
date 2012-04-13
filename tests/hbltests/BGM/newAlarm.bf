USE_MPI_CACHING = 1;
PRINT_DIGITS = -1;

ExecuteAFile ("/Users/apoon/svn/hyphy/HBL/art/BGM/alarm/bayesgraph.ibf");


import_xmlbif("/Users/apoon/svn/hyphy/HBL/art/BGM/alarm/alarm.xml", "Alarm"); 

names = Rows(Alarm);
nodes={};
for (i = 0; i < Abs(Alarm); i=i+1)
{
	nodes[Abs(nodes)] = add_discrete_node (names[i], 2, 0, (Alarm[names[i]])["Levels"]);
}

num_nodes = Abs(nodes);

BGM alarm_bgm = (nodes);

for (rep = 0; rep < 20; rep = rep + 1)
{
	sim = ReadCSVTable ("/Users/apoon/svn/hyphy/HBL/art/BGM/alarm/censor25pc/alarm1000_25pc_"+rep+".csv", 0);
	
	attach_data("alarm_bgm", sim, 250, 0, 100);
	
	result0 = order_MCMC ("alarm_bgm", 100000, 100000, 100);
	write_edgelist("/Users/apoon/svn/hyphy/HBL/art/BGM/alarm/order_"+rep+".edges", result0, num_nodes, 1);
	
	result1 = graph_MCMC ("alarm_bgm", 1000000, 1000000, 100, 0);
	write_edgelist("/Users/apoon/svn/hyphy/HBL/art/BGM/alarm/graph_"+rep+".edges", result1, num_nodes, 1);
}
