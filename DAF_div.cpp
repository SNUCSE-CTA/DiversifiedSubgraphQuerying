#include "DAF.h"

int* marking_cnt;
int cur_k;
int num_k;
class emb_element{
		public:
		int* mapped_vertex;

		emb_element(){
				mapped_vertex = new int[cnt_node_query];
		}
};
vector<emb_element> stored_emb;

bool b_over_time_limit = false;
#ifdef WINDOWS
bool b_search_end = false;
#endif
void* timer(void* args)
{
//	int time_limit_milli = 1000 * 60 * 10; //10 minutes
	int time_limit_milli = 1000 * 10; //10 minutes
	int sleep_time_sec = 1; //1 second
	auto begin = chrono::high_resolution_clock::now();
	
	while(true)
	{
#ifdef WINDOWS
		if (b_search_end)
			return nullptr;
#endif
		auto now = chrono::high_resolution_clock::now();
		chrono::duration<double, milli> elapsed = now - begin;

		if(elapsed.count() > time_limit_milli)
		{
			b_over_time_limit = true;
			return nullptr;
		}
#ifdef WINDOWS
		this_thread::sleep_for(std::chrono::milliseconds(sleep_time_sec * 1000));
#else
		sleep(sleep_time_sec);
#endif
	}

	return nullptr;
}

//for dynamic programming between a directed acyclic graph and a graph
void BFS_DAG();
void print_DAG(); //for debug
void print_CPI(); //for print CPI; TODO: rename function
//dynamic programming between a DAG and a graph
void topDownInitial();
void bottomUpIterate();
void topDownIterate();
void adjacencyListConstruction();
//ancestors for failing set

//===== some key counts ===============================

int * label_degree_to_node;

inline void preprocessDataGraph(string datagraphFile) {
    //preprocess the data graph to extract some key info which is to be used in the next function
	ifstream fin(datagraphFile);
	string line;
	getline(fin, line);
	vector<string> v;
	split(line, v, ' ');
	cnt_node = atoi(v[2].c_str()); //get count_node_number

	nodes_label = new int [cnt_node]; 

	while (getline(fin, line)) {
		if (line.at(0) == 'v') {
			split(line, v, ' ');
			int label = atoi(v[2].c_str());
			if (label > largest_label)
				largest_label = label;
			unique_label_set.insert(label);
		}
		if (line.at(0) == 'e')
			count_edge_number++;
	}

	cnt_unique_label = unique_label_set.size();
	label_freqency = new int [ cnt_unique_label + 1 ];
	memset(label_freqency, 0, sizeof(int) * (cnt_unique_label + 1));
}

inline void readAndProcessDataGraph(string datagraphFile) {

    /* Function readAndProcessDataGraph(string datagraphFile)
    * 1. For vertices, calculate nodes_label[v_id] and label_frequency[actual_label]
    * 2. For edges, calculate the adjacent list nodes[v_id] and data_edge_matrix[v1_id*size+v2_id]
    *    data_edge_matrix[i*size+j] = 1 if there is an edge (i,j). data_edge_matrix[i*size+j] = 0 otherwise
    * 3. For each vertex, store its degree degree_data[v_id] and core_number_data[v_id]
    * 4. For each vertex, create nodes_data[sum_degree], which stores all vertices' adjacent lists, 
    *    nodes_info[v_id], which stores the starting position of v_id's adjacent list in nodes_data
    *    nodes_to_label_info[v_id*(cnt_unique_label+1)+cur_label], which stores the starting positions and number of 
    *    v_id's neighbors with cur_label in nodes_data
    * 5. Build label_degree_to_node[v_id], a sorted array of all vertices by their labels and then ascending order of degree, 
    *    and label_deg_label_pos[label], which stores a pair of 
    *    <the largest degree of this label, the end position (not included) of this label's largest degree >
    * 6. Build the NLF lightweight check
    * 7. Create an |V|-sized array of HashTable only if |V| is larger than the size of edge matrix
    * 8. Compute the maximum neighbor degree MAX_NB_degree_data[v_id] for each vertex
    */

	/* The input file should be in the node and edge format
	 * This function read the data graph and store it in two formats: node, edge and adjacent list
	 * nodes[]: the array store a node's adjacent list
	 * edges[]: store all edges
	 * node_label: store the label of the nodes
	 */

	// used new label to node index with the degree sensitive
	//initialize the array to all 0;

	vector<vector<int> > nodes; //locally in this function, such that it'll be released to save space
	nodes.resize( cnt_node );

	transferredLabel = new int[largest_label + 1];
	memset(transferredLabel, 0, sizeof(int) * (largest_label + 1));
	int label_index = 1; //start from 1, because 0 means not found the label

	nodes_info = new int[cnt_node + 1];
	nodes_to_label_info = new pair<int, int> [cnt_node * (cnt_unique_label + 1)];

	int node_id = 0;

	//===== read data from file ================
	ifstream fin(datagraphFile);
	string line;
	vector<string> v;
	while (getline(fin, line)) {
		if (line.at(0) == 'v') { //this is node
			split(line, v, ' ');
			int label = atoi(v[2].c_str()); // get node label
			int actual_label = transferredLabel[label]; // transfer the label to the actual label we can use
			if (actual_label == 0) { //means this label is first time being transfer
				actual_label = label_index;
				transferredLabel[label] = label_index;
				label_index++; // maintain the label index. notice that 0 means not found here, therefore 0 is not a valid label
			}
			nodes_label[node_id] = actual_label;
			node_id++;
			label_freqency[actual_label] ++;
		}

		if (line.at(0) == 'e') {
			split(line, v, ' ');
			int left_node = atoi(v[1].c_str());
			int right_node = atoi(v[2].c_str());
			if (left_node != right_node){//put the nodes into the adjacent list
				nodes[left_node].push_back(right_node);
				nodes[right_node].push_back(left_node);
			}
		}
	}
	fin.close();
	//===================

	//====== build the degree and core number array ===============
	core_number_data = new int [cnt_node];
	degree_data = new int[cnt_node];
	MAX_DEGREE = 0;
	int sum_degree = 0;
	for (int i = 0; i < cnt_node; i++) {
		int degree = nodes[i].size();
		degree_data[i] = degree;
        //std::cout<<"degree_data["<<i<<"]: "<<degree<<std::endl;
		core_number_data[i] = degree;
		sum_degree += degree;
		if (degree > MAX_DEGREE)
			MAX_DEGREE = degree;
	}
	//==============================================

	nodes_data = new int [sum_degree];

	int cur_array_index = 0;
	int cur_matrix_index = 0;

	for (int i = 0; i < cnt_node; i++) {

		//========================= put this node's sorted adjacent list into the array ================================================
		if (nodes[i].size() == 0){  //deal with isolated nodes, specially for yeast
//			cerr << i <<  " => degree zero node!!" << endl;
			nodes_info[i] = cur_array_index;
			continue;
		}
		nodes_info[i] = cur_array_index; // indicate the starting position of node i's adjacent list in "nodes_data"
		//stable_sort(nodes[i].begin(), nodes[i].end(), sortByLabelAndDegree);//sort by label and then ascending order of degree
		sort(nodes[i].begin(), nodes[i].end(), sortByLabelAndDegree);//sort by label and then ascending order of degree
		copy (nodes[i].begin(), nodes[i].end(), nodes_data + cur_array_index  );
		cur_array_index += nodes[i].size(); // maintain the index
		//=====================================================================================

		int cur_label = nodes_label[ nodes[i][0] ];
		int cur_count = 1;

		if (nodes[i].size() == 1) { //special case: there is only one node in the adjacent list
			nodes_to_label_info[i * (cnt_unique_label + 1) + cur_label] = make_pair(cur_matrix_index, cur_count);
			cur_matrix_index += cur_count;
			continue;
		}

		for (int j = 1; j < nodes[i].size(); j++) {

			int this_label = nodes_label[nodes[i][j]];

			if (this_label == cur_label)
				cur_count++;
			else {
				nodes_to_label_info[i * (cnt_unique_label + 1) + cur_label] = make_pair(cur_matrix_index, cur_count);
				cur_matrix_index += cur_count;
				cur_label = this_label;
				cur_count = 1;
			}

		}

		nodes_to_label_info[i * (cnt_unique_label + 1) + cur_label] = make_pair(cur_matrix_index, cur_count);
		cur_matrix_index += cur_count;
		//=====================================================================================

	} //end for

	nodes_info[cnt_node] = sum_degree; //the end position for the last node of the nodes_info

	//=========== build a structure for the "label to nodes with degree" array =============
	int last_label = 0;
	int last_degree = 0;
	label_degree_to_node = new int[cnt_node];
	for (int i = 0; i < cnt_node; i++)
		label_degree_to_node[i] = i;

	//stable_sort(label_degree_to_node, label_degree_to_node + cnt_node, sortByLabelandDegree);
	sort(label_degree_to_node, label_degree_to_node + cnt_node, sortByLabelandDegree);

	//<the largest degree of this label, the end position (not included) of this label's largest degree >
	label_deg_label_pos = new pair<int, int> [ cnt_unique_label + 1 ];

	degree_array.resize(cnt_node);

	for (int i = 0; i < cnt_node; i++)
		degree_array [i] = degree_data [  label_degree_to_node [i] ];

	for (int i = 0; i < cnt_node; i++) {
		int v = label_degree_to_node[i];
		int label = nodes_label[ v ];
		int degree = degree_data[ v ];
		if (i == 0) {
			label_deg_label_pos[ 0 ] = make_pair(0, 0); //actully, "0" is not used as label, so not necessary
			last_label = label;
			last_degree = degree;
		} else {
			if (label != last_label) //deal with a new label,
				label_deg_label_pos[ last_label ] = make_pair(last_degree, i);
			last_label = label;
			last_degree = degree;
		}
	}

	label_deg_label_pos[ last_label ] = make_pair(last_degree, cnt_node);//degree value is meaningless here.
	//========================================================================================

	NEC_mapping = new int[(cnt_unique_label + 1) * MAX_QUERY_NODE]; //label and the parent node id
	memset(NEC_mapping, 0, sizeof(int) * (cnt_unique_label + 1) * MAX_QUERY_NODE);

	NEC_mapping_pair = new NEC_element[(cnt_unique_label + 1) * MAX_QUERY_NODE]; // <label, parent node id>
	NEC_set_array = new NEC_set_array_element[(cnt_unique_label + 1) * MAX_QUERY_NODE]; // <parent_id, sum>

	//================= build the NFL lightweight check =====================
	NLF_size = (cnt_unique_label + 1) / SIZEOF_INT + 1;
	NLF_array = new int [ NLF_size ]; // the array for query graph
	NLF_check = new int [cnt_node * NLF_size];
	memset(NLF_check, 0, sizeof(int) * NLF_size * cnt_node);
	for (int i = 0; i < nodes.size(); i++)
		for (int j = 0; j < nodes[i].size(); j++){
			int label = nodes_label[ nodes[i][j] ];
			int idx = NLF_size - 1 - label / SIZEOF_INT;
			int pos = label % SIZEOF_INT;
			NLF_check[i  * NLF_size + idx] |= (1 << pos);
		}
	//=====================================================================================
	//============== initialize the MAX neighbor degree of data node =====================
	MAX_NB_degree_data = new int [cnt_node];
	for (int i = 0; i < cnt_node; i++){
		int max_degree = 0;
		for (int j = 0; j < nodes[i].size(); j++){
			int node = nodes[i][j];
			if (degree_data[node] > max_degree)
				max_degree = degree_data[node];
		}
		MAX_NB_degree_data[i] = max_degree;
	}
	//=====================================================================================
}

inline void dataGraphAnalysis(string datagraphFile) {

	/* The input file should be in the node and edge format
	 * This function read the data graph and store it in two formats: node, edge and adjacent list
	 * nodes[]: the array store a node's adjacent list
	 * edges[]: store all edges
	 * node_label: store the label of the nodes
	 */
    int AVERAGE_DEGREE_DATA_SUM = 0;

	{
		ifstream fin(datagraphFile);
		string line;
		getline(fin, line);
		vector<string> v;
		split(line, v, ' ');

		cnt_node = atoi(v[2].c_str()); //get count_node_number

		while (getline(fin, line)) {

			if (line.at(0) == 'v') {
				split(line, v, ' ');
				int label = atoi(v[2].c_str());
				if (label > largest_label)
					largest_label = label;
				unique_label_set.insert(label);
			}

			if (line.at(0) == 'e')
				count_edge_number++;
		}

		cnt_unique_label = unique_label_set.size();

		label_freqency = new int [ cnt_unique_label + 1 ];
		memset(label_freqency, 0, sizeof(int) * (cnt_unique_label + 1));

	}

	vector<vector<int> > nodes;
	nodes.resize( cnt_node );
	// used new label to node index with the degree sensitive
	//initilize the array to all 0;
	transferredLabel = new int[largest_label + 1];
	memset(transferredLabel, 0, sizeof(int) * (largest_label + 1));
	int label_index = 1; //start from 1, because 0 means not found the label
	ifstream fin(datagraphFile);
	string line;
	nodes_info = new int[cnt_node + 1];
	nodes_to_label_info = new pair<int, int> [cnt_node * (cnt_unique_label + 1)];

	int node_id = 0;
	vector<string> v;

	while (getline(fin, line)) {

		if (line.at(0) == 'v') { //this is node
			split(line, v, ' ');
			int label = atoi(v[2].c_str()); // get node label
			int actual_label = transferredLabel[label]; // transfer the label to the actual label we can use

			if (actual_label == 0) { //means this label is first time being transfer
				actual_label = label_index;
				transferredLabel[label] = label_index;
				label_index++; // maintain the label index. notice that 0 means not found here, therefore 0 is not a valid label
			}

			//nodes_label.push_back(actual_label); // now that we use the transferred actual label
			nodes_label[node_id] = actual_label;
			node_id++;
			label_freqency[actual_label] ++;
		}

		if (line.at(0) == 'e') {
			split(line, v, ' ');
			int left_node = atoi(v[1].c_str());
			int right_node = atoi(v[2].c_str());
			if (left_node != right_node){	//put the nodes into the adjacent list
				nodes[left_node].push_back(right_node);
				nodes[right_node].push_back(left_node);
			}
		}
	}

	//====== build the degree vector ===============
	core_number_data = new int [cnt_node];
//	degree_data.resize(cnt_node);
	degree_data = new int[cnt_node];
	MAX_DEGREE = 0;
	int node_degree_largest;
	for (int i = 0; i < cnt_node; i++) {
		degree_data[i] = nodes[i].size();
		core_number_data[i] = nodes[i].size();
		if (degree_data[i] > MAX_DEGREE) {
			MAX_DEGREE = degree_data[i];
			node_degree_largest = i;
		}
	}
	//==============================================


	int sum_degree_one = 0;
	int sum_degree_zero = 0;

	int sum_degree = 0;
	for (int i = 0; i < cnt_node; i++) {
		int degree = nodes[i].size();
		if (degree == 1)
			sum_degree_one ++;
		if (degree == 0)
			sum_degree_zero ++;
		sum_degree += degree;
	}

	AVERAGE_DEGREE_DATA_SUM += (double) sum_degree / cnt_node;

	cout << "Data graph info ******* " << datagraphFile << " ********" << endl;
	cout << "Data graph nodes: " << cnt_node << endl;
	cout << "Data graph edges: " << count_edge_number << endl;
	cout << "Largest label is " << largest_label << endl;
	cout << "Unique labels => " << cnt_unique_label << endl;
	cout << "Average degree is " << AVERAGE_DEGREE_DATA_SUM << endl;
	cout << "Largest degree is " << MAX_DEGREE << endl;
	cout << "Leaf node number: " << sum_degree_one << endl;
	cout << "Isolated nodes  : " << sum_degree_zero << endl << endl;
}

inline void coreDecomposition_data(){
	//core-decomposition for the data graph
	//begin starting the core-decomposition, core number is the degree number

	int* bin = new int[MAX_DEGREE + 1];
	memset(bin, 0, sizeof(int) * ( MAX_DEGREE + 1) );

	for (int i = 0; i < cnt_node; i ++)
		bin[ core_number_data [i] ] ++;

	int start = 0;
	int num;

	for (int d = 0; d <= MAX_DEGREE; d++){
		num = bin[d];
		bin[d] = start;
		start += num;
	}

	int* pos = new int [cnt_node];
	int* vert = new int [cnt_node];

	for (int i = 0; i < cnt_node; i++){
		pos[i] = bin[ core_number_data[i] ];
		vert[ pos[i] ] = i;
		bin[ core_number_data[i] ] ++;
	}

	for (int d = MAX_DEGREE; d > 0; d --)
		bin[d] = bin[d-1];
	bin[0] = 0;

	for (int i = 0; i < cnt_node; i++){

		int v = vert[i];

		for (int j = nodes_info[v]; j < nodes_info[v + 1]; j++){

			int u = nodes_data[j];

			if (core_number_data[u] > core_number_data[v]){

				int du = core_number_data[u];
				int pu = pos[u];

				int pw = bin[du];
				int w = vert[pw];

				if (u != w){	//if not the same node, switch the position of the two nodes.
					pos[u] = pw;
					pos[w] = pu;
					vert[pu] = w;
					vert[pw] = u;
				}

				bin[du] ++;
				core_number_data[u]--;
			}
		}
	}

	delete[] bin;
	delete[] pos;
	delete[] vert;
}

//=========================================================================================

inline void readQueryGraph() {

	int m;
	int matrix_index = 0;
	MAX_DEGREE_QUERY = 0;

	//resize only if required
	if (sum_degree_cur > MAX_sum_degree_cur){
		MAX_sum_degree_cur = sum_degree_cur;
		query_nodes_array.resize(sum_degree_cur * 2);
		core_tree_node_child_array.resize(sum_degree_cur);
		core_tree_node_nte_array.resize(sum_degree_cur);
	}

	for (int i = 0; i < cnt_node_query; i++){

		query_nodes_array_info[i] = matrix_index;

		fin_query >> m >> label_cur >> degree_cur;

		if (degree_cur > MAX_DEGREE_QUERY)
			MAX_DEGREE_QUERY = degree_cur;

		nodes_label_query[i] = transferredLabel[label_cur];
		node_degree_query[i] = degree_cur;
		core_number_query[i] = degree_cur;

		for (int j = 0; j < degree_cur; j ++){
			fin_query >> query_nodes_array[matrix_index];
			query_neighbor_index[i][query_nodes_array[matrix_index]] = j;
			matrix_index ++;
		}
	}
	query_nodes_array_info[cnt_node_query] = matrix_index;
}

inline void query_analysis() {

	cout << "========= Analyzing query set: " << querygraphFileFolder << " ==============="<< endl;

	fin_query.open(querygraphFileFolder);

	char c;
	int m;

	double sum_average_degree = 0;
	double sum_max_degree = 0;

	double max_average_degree = 0;
	double max_max_degree = 0;

	double sum_node_number = 0;

	for (int i = 0; i < count_query_file; i++) {

		fin_query >> c >> m >> cnt_node_query >> sum_degree_cur;

		double average_degree = (double) sum_degree_cur / (double) cnt_node_query;

		sum_node_number += cnt_node_query;
		sum_average_degree += average_degree;

		if (max_average_degree < average_degree)
			max_average_degree = average_degree;

		double max_degree = 0;

		for (int j = 0; j < cnt_node_query; j++) {
			fin_query >> m >> label_cur >> degree_cur;
			if (degree_cur > max_degree)
				max_degree = degree_cur;

			for (int x = 0; x < degree_cur; x ++)
						fin_query >> m;
		}

		sum_max_degree += max_degree;

		if (max_max_degree < max_degree)
			max_max_degree = max_degree;
	}

	fin_query.close();

	cout << "Average number of nodes per query is " << sum_node_number / (double)count_query_file << endl;
	cout << "Average degree of all query is " <<  sum_average_degree /  (double)count_query_file << endl;
	cout << "Average max degree of all query is " << sum_max_degree  / (double) count_query_file << endl;
	cout << "Max average degree of all query is " << max_average_degree << endl;
	cout << "Max max degree of all query is " << max_max_degree  << endl;
	cout << endl;
}

inline void coreDecomposition_query(){

	//begin starting the core-decomposition, core number is the degree number
	int * bin = bin_query;	//	int bin [MAX_DEGREE_QUERY + 1];
	int * pos = pos_query;	//	int pos [cnt_node_query];
	int * vert = vert_query;//	int vert [cnt_node_query];

	memset(bin, 0, sizeof(int) * ( MAX_DEGREE_QUERY + 1) );

	for (int i = 0; i < cnt_node_query; i ++)
		bin[ core_number_query [i] ] ++;

	int start = 0;
	int num;

	for (int d = 0; d <= MAX_DEGREE_QUERY; d++){
		num = bin[d];
		bin[d] = start;
		start += num;
	}


	for (int i = 0; i < cnt_node_query; i++){
		pos[i] = bin[ core_number_query[i] ];
		vert[ pos[i] ] = i;
		bin[ core_number_query[i] ] ++;
	}

	for (int d = MAX_DEGREE_QUERY; d > 0; d --)
		bin[d] = bin[d-1];
	bin[0] = 0;

	for (int i = 0; i < cnt_node_query; i++){

		int v = vert[i];

		for (int j = query_nodes_array_info[v]; j < query_nodes_array_info[v + 1]; j ++){

			int u = query_nodes_array[j]; // nodes_query[v][j];

			if (core_number_query[u] > core_number_query[v]){

				int du = core_number_query[u];
				int pu = pos[u];

				int pw = bin[du];
				int w = vert[pw];

				if (u != w){	//if not the same node, switch the position of the two nodes.
					pos[u] = pw;
					pos[w] = pu;
					vert[pu] = w;
					vert[pw] = u;
				}

				bin[du] ++;
				core_number_query[u]--;
			}
		}
	}

}

//======================================FUNCTIONS EXCLUSIVELY FOR NORMAL INPUT ====================================================

inline void constructNECMapping(int a_node)
{
	//for all of node i's children
	for (int j = query_nodes_array_info[a_node]; j < query_nodes_array_info[a_node + 1]; j ++){

	    int child = query_nodes_array[j];

		if (core_number_query[child] < 2){ // the child node is not in the 2-core

		//two cases here, the NEC node or a residual tree

		    if (node_degree_query[child] == 1){ //degree is one ==> NEC node
		        //============ CASE ONE: ONE-DEGREE NODES => NEC nodes =====================

			    int label = nodes_label_query[child];

				if (NEC_mapping[label * MAX_QUERY_NODE + a_node] == 0) {
					NEC_mapping_pair[NEC_mapping_pair_index ++] = NEC_element(label, a_node, child);// child is the representative node
					NEC_map [child] = child;//NEC map
					NEC_mapping_Actual[label * MAX_QUERY_NODE + a_node] = child;
				} else {
					NEC_map [child] = NEC_mapping_Actual[label * MAX_QUERY_NODE + a_node];//NEC map
				}

				NEC_mapping[label * MAX_QUERY_NODE + a_node] ++; // the label with parent being i, nec_count ++

			} else {
				//============ CASE TWO: NORMAL CASE, THE QUERY TREE ================
				// extract the query tree for extra region candidate extraction, based on DFS
				// also give a DFS-based query sequence at the same time

				int * dfs_stack = dfs_stack_query;
				int dfs_stack_index = 0;

				visited_for_query[a_node] = 1; // this is the start node's parent node (a marked node)
				visited_for_query[child] = 1; // this is the start node

				dfs_stack[dfs_stack_index ++] = child;
				residual_tree_match_seq[residual_tree_match_seq_index ++] = child;

				tree_node_parent[child] = a_node;

				while (dfs_stack_index != 0) {

					int current_node = dfs_stack[dfs_stack_index - 1];
					dfs_stack_index--;

					int added_child = 0;

					for (int m = query_nodes_array_info[current_node]; m < query_nodes_array_info[current_node + 1]; m ++){

						int child_node = query_nodes_array[m];

						if (!visited_for_query[child_node]) {

							visited_for_query[child_node] = 1;

							//======== special treatment here: if a node is a leaf (degree being 1), then put it into nec node set
							if (node_degree_query[child_node] == 1){

								int label = nodes_label_query[child_node];

								if (NEC_mapping[label * MAX_QUERY_NODE + current_node] == 0) {
									NEC_mapping_pair[NEC_mapping_pair_index ++] = NEC_element(label, current_node, child_node);// child is the repesentive node
									NEC_map [child_node] = child_node;//NEC map
									NEC_mapping_Actual[label * MAX_QUERY_NODE + current_node] = child_node;
								}
								else{
									NEC_map [child_node] = NEC_mapping_Actual[label * MAX_QUERY_NODE + current_node];//NEC map
								}
								NEC_mapping[label * MAX_QUERY_NODE + current_node]++; // the label with parent being i, nec_count ++
								continue;
							}
							//===========================================================
							tree_node_parent[child_node] = current_node;
							added_child ++;
							dfs_stack[dfs_stack_index ++] = child_node;
							residual_tree_match_seq[residual_tree_match_seq_index ++] = child_node;
						}

						if (added_child == 0)//this information is recorded for extracting the matching sequence for the tree matching sequence.
							residual_tree_leaf_node[residual_tree_leaf_node_index ++] = make_pair(current_node, 0);
					}
				}
			}
		}
	}
}

inline void extractResidualStructures(){

	residual_tree_match_seq_index = 0;
	residual_tree_leaf_node_index = 0;
	NEC_mapping_pair_index = 0;

	memset(NEC_map, -1, sizeof(int) * cnt_node_query );

	memset(visited_for_query, 0, sizeof(char) * cnt_node_query);

    if(isTree)
    {
        //exit(1);
        //constructNECMapping(root_node_id);
    }
    else
    {
	    for (int i = 0; i < cnt_node_query; i++) {//for each node in the query

		    if (core_number_query[i] < 2)//not in the two-core
			    continue;

		    //now i must be a 2-core node => next, we check whether i is a articulation node
            //constructNECMapping(i);
	    }
    }

	//================ construct the NEC set by label: each label is with a vector which contains many NECs with this label.=========
	sort(NEC_mapping_pair, NEC_mapping_pair + NEC_mapping_pair_index, sort_by_NEC_label);
	int last_label;
	NEC_set_index = 0;
	NEC_set_by_label_index.clear();
	int sum;
	if (NEC_mapping_pair_index == 1){
		NEC_element & nec_ele = NEC_mapping_pair[0];
		int label = nec_ele.label;
		int parent_id = nec_ele.parent_id;
		int represent_child = nec_ele.represent_node;
		sum = NEC_mapping[label * MAX_QUERY_NODE + parent_id];
		NEC_mapping[label * MAX_QUERY_NODE + parent_id] = 0; //reset it back to 0
		NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
		NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
		NEC_set_by_label_index.push_back(make_pair(-1, NEC_mapping_pair_index)); // redundant element to set the end
	} else {
		for (int i = 0; i < NEC_mapping_pair_index; i++) {

			NEC_element & nec_ele = NEC_mapping_pair[i];

			int label = nec_ele.label;
			int parent_id = nec_ele.parent_id;
			int represent_child = nec_ele.represent_node;
			sum = NEC_mapping[label * MAX_QUERY_NODE + parent_id];
			NEC_mapping[label * MAX_QUERY_NODE + parent_id] = 0; //reset it back to 0

			if (i == 0) {
				NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
				NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
				last_label = label;
				continue;
			} else if (i == NEC_mapping_pair_index - 1) {
				if (label != last_label)
					NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
				NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
				NEC_set_by_label_index.push_back(make_pair(-1, NEC_mapping_pair_index)); // redunant element to set the end
				continue;
			}

			if (label != last_label) {
				NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
				last_label = label;
			}

			NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
		}
	}

#ifdef	OUTPUT_EXTRA_INFO

	int sum_node = 0;

	if (NEC_mapping_pair_index != 0){
		for (int i = 0; i < NEC_set_by_label_index.size() - 1; i++) {
			int label = NEC_set_by_label_index[i].first;
			int start = NEC_set_by_label_index[i].second;
			int end = NEC_set_by_label_index[i + 1].second;

			for (int j = start; j < end; j++) {
				int parent_id = NEC_set_array[j].parent_id;
				int sum = NEC_set_array[j].sum;
				sum_node += sum;
				cerr << "label :" << label << " => parent id " << parent_id << " \t sum => " << sum
						<< "\t representative node is " << NEC_set_array[j].represent_node<< endl;
			}
		}
	}

	cerr << "NEC classes contained: " << NEC_mapping_pair_index << " classes with " << sum_node << " nodes " << endl;
	cerr << "Query trees with sum node: " << residual_tree_match_seq_index
				<< " and tree leaf index is " << residual_tree_leaf_node_index << endl;
    if(isTree)
    {
	    cerr << "Nodes in tree: ";
	    for (int i = 0; i < residual_tree_match_seq_index; i++){
		    cerr << residual_tree_match_seq[i] << " ";
	    }
	    cerr << endl;
    }
#endif
}

inline int start_node_selection(){

	double least_ranking = DBL_MAX;
	int start_node = -1;
	double ranking;
	int label;
	int degree;

	for (int i = 0; i < cnt_node_query; i++){
		if (core_number_query[i] < 2 && isTree == false)	//root node must be selected from the core structure
			continue;

		label = nodes_label_query[i];
		degree = node_degree_query[i];

		//binary search used here
		int s = label_deg_label_pos[ label - 1 ].second;
		int end = label_deg_label_pos[ label ].second;
		vector<int>::iterator pos = lower_bound( degree_array.begin() + s , degree_array.begin() + end, degree);
		int start = pos - degree_array.begin();

		ranking = (double)(end - start) /(double)degree ;

		if (ranking < least_ranking){
			least_ranking = ranking;
			start_node = i;
		}
	}
	return start_node;
}

inline void BFS_DAG() {

	/*
	 * output : true_leaf_nodes, simulation_sequence_array, level_to_sequence which maps a level to a segment in the sequence
	 */

	core_tree_node_child_array_index = 0;
	core_tree_node_nte_array_index = 0;

    bfs_sequence_index = 0; // 20170503
	exploreCRSequence_indx = 0;

	char * popped = visited_for_query;			//popped[i] = 1 if i-node have popped from queue.
	memset(popped, 0, sizeof(char) * cnt_node_query);

	int * visited = visited_int_for_query;			//visited[i] = level if i-node had pushed into queue, where level is it's BFS-level.
	memset(visited, 0, sizeof(int) * cnt_node_query);

	int * queue_array = queue_array_query;

	queue_array[0] = root_node_id;

	int pointer_this_start = 0;
	int pointer_this_end = 1;
	int pointer_next_end = 1;
	int current_level = 1; //initial level starts at 1

	simulation_sequence_index = 0;
	level_index.clear();
	visited[root_node_id] = 1;
	BFS_level_query[root_node_id] = 1;
	BFS_parent_query[root_node_id] = -1;

    if(use_failing_set)
    {
        if(cnt_node_query != prev_cnt_node_query){
            zero_vector = bit_vector(cnt_node_query, 0);
            one_vector  = bit_vector(cnt_node_query, 1);
            prev_cnt_node_query = cnt_node_query;
        }
    }

	for(int i = 0; i < cnt_node_query; i++)
	{
		DAG_child_query[i] = new int[node_degree_query[i]];
		memset(DAG_child_query[i], -1, sizeof(int) * node_degree_query[i] );
		DAG_parent_query[i] = new int[node_degree_query[i]];
		memset(DAG_parent_query[i], -1, sizeof(int) * node_degree_query[i] );
		DAG_child_query_parent_index[i] = new int[node_degree_query[i]];
		memset(DAG_child_query_parent_index[i], -1, sizeof(int) * node_degree_query[i] );

		//DAG_ancestor[i] = bit_vector(cnt_node_query, 0);
	}
	memset(DAG_child_query_size, 0, sizeof(int) * cnt_node_query);
	memset(DAG_parent_query_size, 0, sizeof(int) * cnt_node_query);
    //<sort by label_freq first, and sort by degree for each label group>
	//sort label by label frequency
	label_frequency_rank = new int [ cnt_unique_label + 1 ]; //later, these lines should be moved to preprocessData...()
	int* temp_sort_array = new int [ cnt_unique_label + 1 ]; 
	for(int i = 1; i < cnt_unique_label+1; i++) temp_sort_array[i] = i; //label uses [1-cnt_unique_label+1]
	//label_frequency_rank[i] contains label i's rank in terms of the label_frequency. high rank has high frequency; i.e., label_frequency_rank[0] has the biggest label_frequency.
	sort(temp_sort_array+1, temp_sort_array + cnt_unique_label+1, sortByLabelFrequency_Label_inc);
	for(int i =1; i < cnt_unique_label+1; i++) label_frequency_rank[ temp_sort_array[i] ] = i;
    //<\sort by label_freq first, and sort by degree for each label group>

	while (true) {

		int start = simulation_sequence_index;
        //<sort by label_freq first, and sort by degree for each label group>
		sort(queue_array+pointer_this_start, queue_array+pointer_this_end, sortByDegree_Query_dec);
		//by sorting the array using label_frequency_rank, where distinct label has distinct rank, and then only the array can be grouped by labels.
		stable_sort(queue_array+pointer_this_start, queue_array+pointer_this_end, sortByLabelFrequencyRank_Query_inc); 
        //<\sort by label_freq first, and sort by degree for each label group>
		while (pointer_this_start != pointer_this_end) { // queue not empty

			int current_node = queue_array[pointer_this_start];
			pointer_this_start++;
			popped[current_node] = 1;

			int start = query_nodes_array_info[current_node];
			int end = query_nodes_array_info[current_node + 1];
			//cout<<"BFS_DAG. current_node: "<<current_node<<", (start, end): ("<<start<<", "<<end<<")"<<endl;	

			for (int i = start; i < end; i ++){

				int childNode = query_nodes_array[i];
                //cout<<"childNode "<<i<<": "<<childNode<<endl;
				if (popped[childNode] == 0) //childNode is not current_node's parent
				{
					DAG_child_query[current_node][ DAG_child_query_size[current_node] ] = childNode;
					DAG_child_query_parent_index[current_node][ DAG_child_query_size[current_node] ] = DAG_parent_query_size[childNode];

					DAG_parent_query[childNode][ DAG_parent_query_size[childNode] ] = current_node;

					//DAG_ancestor[childNode][current_node] = 1;	       //base case
					//DAG_ancestor[childNode] |= DAG_ancestor[current_node]; //propagate

					DAG_child_query_size[current_node]++;
					DAG_parent_query_size[childNode]++;
				}

                if (visited[childNode] == 0) //this child node has not been visited.
                {

					visited[childNode] = current_level + 1; //parent node's level plus one

					queue_array[pointer_next_end] = childNode;
					pointer_next_end++;

					BFS_level_query[childNode] = current_level + 1;
					BFS_parent_query[childNode] = current_node;

					if (core_number_query[childNode] < 2)
						continue;
				}
			}

			simulation_sequence[simulation_sequence_index] = current_node;
			simulation_sequence_index ++;
		}

		int end = simulation_sequence_index;

		level_index.push_back(make_pair(start, end));
        //cout <<"BFS_DAG. current_level: "<<current_level<<", start: "<<start<<", end: "<<end<<endl;
		for (int i = start; i <= end-1; i++){ //for (int i = end - 1; i >= start; i--){
			int node = simulation_sequence[i];
            //if(node == root_node_id or NEC_map[node] == -1){ //root_nodeid can be NEC node!
            if( NEC_map[node] == -1){ //root_nodeid can be NEC node!
                bfs_sequence_index++;
            }
			if (core_number_query[node] < 2)
				continue;
			exploreCRSequence_indx++;
		}

		if (pointer_next_end == pointer_this_end) //no node has been pushed in
			break;

		pointer_this_start = pointer_this_end;
		pointer_this_end = pointer_next_end;

		current_level++;

	}
}

//<bad_array_new_length>
inline pair<int, int> nodesToLabelInfo(int aVertex, int aLabel)
{
    int low = nodes_info[aVertex];
    int high = nodes_info[aVertex + 1];

    int start = -1;
    int count = 0;
    //linear search to the start and end
    for (int i = low; i < high; i++)
    {
        if (nodes_label[nodes_data[i]] == aLabel)
        {
            if (start == -1)
                start = i;

            count++;
        }
    }
    pair<int, int> ret(start, count);
    return ret;
}
//<\bad_array_new_length>

inline void topDownInitial () {
    //cout << "root node: " << root_node_id << endl;
    memset(flag_prelin_char, 0, sizeof(char)*cnt_node);
	//1. initialize candidates of root vertex
	{
		NodeIndexUnit & root_node_unit = indexSet[root_node_id];
		int label = nodes_label_query[root_node_id];
		int degree = node_degree_query[root_node_id];

		//1.1. make NLF, core_degree filter, and max_nb_degree filter
		int max_nb_degree = 0;
		//core_degree
		int core = core_number_query[root_node_id];
        //cout << "core: " << core << endl;
		//NLF
		int first = query_nodes_array_info[root_node_id];
		memset(NLF_array, 0, sizeof(int) * NLF_size);
		for (int j = first; j < first + degree; j++) {

			int local_label = nodes_label_query[  query_nodes_array[j] ];
			int idx = NLF_size - 1 - local_label / SIZEOF_INT;
			int pos = local_label % SIZEOF_INT;
			NLF_array[idx] |= (1 << pos);

			//max_nb_degree
			int nb_degree = node_degree_query [ query_nodes_array[j] ];
			if ( nb_degree > max_nb_degree) //find the max neighbor degree
				max_nb_degree = nb_degree;
		}
		//end 1.1.

		//1.2. for each v in G with label and passing filters
		int s = label_deg_label_pos[ label - 1 ].second;
		int end = label_deg_label_pos[ label ].second;
		vector<int>::iterator pos = lower_bound( degree_array.begin() + s , degree_array.begin() + end, degree);
		int start = pos - degree_array.begin();

		count_global_temp_array_1 = 0;
        //cout << "can_id: ";
		//for each v in G with label
		for (int j = start; j < end; j++) {

			int can_id = label_degree_to_node[j]; //v
            //cout << can_id << " ";
			if (core_number_data[can_id] < core || max_nb_degree > MAX_NB_degree_data[can_id])
				continue;

			char flag_add = 1;
			for (int pos = NLF_size - 1; pos >= 0; pos--){
				if (NLF_check[(long long)can_id * (long long)NLF_size + (long long)pos] != ( NLF_array[pos] | NLF_check[(long long)can_id * (long long)NLF_size + (long long)pos] )){
					flag_add = 0;
					break;
				}
			}

			if (flag_add)
				root_node_unit.candidates[count_global_temp_array_1 ++] = can_id;

		} // end for
        //cout << endl;
		root_node_unit.size = count_global_temp_array_1;
	}
	//end 1.

	//2. for each query vertex u of q in a top-down fashion
	for (int i = 1; i < simulation_sequence_index; i ++){ // 'i=0' is root, which is already processed above.

		int current_node = simulation_sequence[i]; //u
		int label_cur = nodes_label_query[current_node];
		int degree_cur = node_degree_query[current_node];
		char check_value = 0;	//Cnt = 0

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node)
			continue;

		//2.1. make NLF, core_degree filter, and max_nb_degree filter
		int max_nb_degree = 0;
		//core_degree
		int core_cur = core_number_query[current_node];
		//NLF
		int first = query_nodes_array_info[current_node];
		memset(NLF_array, 0, sizeof(int) * NLF_size);
		for (int j = first; j < first + degree_cur; j++) {
			int local_label = nodes_label_query[  query_nodes_array[j] ];
			int idx = NLF_size - 1 - local_label / SIZEOF_INT;
			int pos = local_label % SIZEOF_INT;
			NLF_array[idx] |= (1 << pos);

			//max_nb_degree
			int nb_degree = node_degree_query [ query_nodes_array[j] ];
			if ( nb_degree > max_nb_degree) //find the max neighbor degree
				max_nb_degree = nb_degree;
		}
		//end 2.1.
		
		//2.2. for each parent p of u in dag
		for(int dag_parent_index = 0; dag_parent_index < DAG_parent_query_size[current_node]; dag_parent_index++){ 

			int parent = DAG_parent_query[current_node][dag_parent_index]; //p
			NodeIndexUnit & parent_node_unit = indexSet[parent];

			//2.2.1. for each candidate vertex v' in p.C
			for (int y = 0; y < parent_node_unit.size; y++) {
				int parent_cand = parent_node_unit.candidates[y]; //v'

				//2.2.1.1. for each vertex v adjacent to v' with label
				pair<int, int> query_result = nodes_to_label_info[parent_cand * (cnt_unique_label + 1) + label_cur];
                //pair<int, int> query_result = nodesToLabelInfo(parent_cand, label_cur);
				for (int z = query_result.first; z < query_result.first + query_result.second; z++) {
					int can_id = nodes_data[z];
					if (flag_prelin_char[can_id] == check_value){ //CHECK: char type can only store 8-bits. overflow check
						flag_prelin_char[can_id] ++; //v.cnt++;
						if (check_value == 0)
							array_to_clean[ to_clean_index ++ ] = can_id;
					}
				}
			}
			check_value ++; //update the check value by one
		}
		//end 2.2.

		//2.3. for each vertex v with v.cnt = Cnt, if v passes filters, add candidates
		NodeIndexUnit & cur_node_unit = indexSet[current_node];

		int cand_index = 0;
		for(int check_index = 0; check_index < to_clean_index; check_index++){ 

			int can_id = array_to_clean[check_index]; //v
			if( flag_prelin_char[can_id] == check_value ){

				//check degree, core, and max_neighbor degree together
				if (degree_data[can_id] < degree_cur || core_number_data[can_id] < core_cur || MAX_NB_degree_data[can_id] < max_nb_degree)
					continue;

				//check lightweight NLF
				char flag_add = 1;
				for (int pos = NLF_size - 1; pos >= 0; pos--){
					if (NLF_check[(long long)can_id * (long long)NLF_size + (long long)pos] != ( NLF_array[pos] | NLF_check[(long long)can_id * (long long)NLF_size + (long long)pos] )){
						flag_add = 0;
						break;
					}
				}

				// lightweight NLF OK
				if (flag_add){
					cur_node_unit.candidates[cand_index] = can_id;
					cand_index++;
				}
			}

		}
		cur_node_unit.size = cand_index;
		//end 2.3.

		//2.4. reset v.cnt for all v in G such that v.cnt > 0	
		while(to_clean_index != 0) //this can erase where to_clean_index==0 because below line pre-decrease to_clean_index
			flag_prelin_char [ array_to_clean[ --to_clean_index ] ] = 0; 
		//end 2.4.

	}//for simulation sequence
	//end 2.
}

inline void bottomUpIterate (){

	//1. for each query vertex u of dag(q) in a bottom-up fashion
	for (int i = simulation_sequence_index-1; i >= 0; i--){ 

		int current_node = simulation_sequence[i]; //u
		int label_cur = nodes_label_query[current_node];

		char check_value = 0;	//Cnt = 0

		if (DAG_child_query_size[current_node] == 0) //No child, No prunning.
			continue;
		
		//1.1. for each child uc of u in dag(q)
		for(int dag_child_index = 0; dag_child_index < DAG_child_query_size[current_node]; dag_child_index++){ 

			int child = DAG_child_query[current_node][dag_child_index]; //uc
			NodeIndexUnit & child_node_unit = indexSet[child];

			//NEC boost
			if (NEC_map[child] != -1 && NEC_map[child] != child) //process only represent node when the child is in NEC
				continue;

			//1.1.1. for each candidate vc in uc.C
			for (int y = 0; y < child_node_unit.size; y++) {
				int child_cand = child_node_unit.candidates[y]; //vc

				//1.1.1.1 for each vertex v adjacent to vc
				pair<int, int> query_result = nodes_to_label_info[child_cand * (cnt_unique_label + 1) + label_cur];
//                pair<int, int> query_result = nodesToLabelInfo(child_cand, label_cur);
				for (int z = query_result.first; z < query_result.first + query_result.second; z++) {
					int can_id = nodes_data[z]; //v
					if (flag_prelin_char[can_id] == check_value){
						flag_prelin_char[can_id] ++;
						if (check_value == 0)
							array_to_clean[ to_clean_index ++ ] = can_id;
					}
				}
			}
			check_value ++; //update the check value by one
		}
		//end 1.1.

		//1.2. for each vertex v in u.C	
		NodeIndexUnit & cur_node_unit = indexSet[current_node];

		for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++){ 
			int can_id = cur_node_unit.candidates[cand_index]; //v
			if( flag_prelin_char[can_id] != check_value ){ //if v.cnt != Cnt then Remove v from u.C
				cur_node_unit.candidates[cand_index] = cur_node_unit.candidates[ --cur_node_unit.size];
				--cand_index;
			}
		}
		//end 1.2.

		//1.3. reset v.cnt = 0 for all v in G such that v.cnt > 0
		while(to_clean_index != 0) //reset v.cnt = 0 for all v in G s.t. v.cnt > 0;
			flag_prelin_char [ array_to_clean[ --to_clean_index ] ] = 0;
		//end 1.3.

	}//for simulation sequence
}
	
inline void topDownIterate (){

	//1. for each query vertex u of dag(q) in a top-down fashion
	for (int i = 1; i < simulation_sequence_index; i++){ //root (0) has no parent

		int current_node = simulation_sequence[i]; //u
		int label_cur = nodes_label_query[current_node];

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node)
			continue;

		char check_value = 0;	//Cnt = 0

		//1.1. for each parent up of u in dag(q)
		for(int dag_parent_index = 0; dag_parent_index < DAG_parent_query_size[current_node]; dag_parent_index++){ 

			int parent = DAG_parent_query[current_node][dag_parent_index]; //up
			NodeIndexUnit & parent_node_unit = indexSet[parent];

			//1.1.1. for each candidate vp in up.C
			for (int y = 0; y < parent_node_unit.size; y++) {
				int parent_cand = parent_node_unit.candidates[y]; //vp

				//1.1.1.1 for each vertex v adjacent to vp
				pair<int, int> query_result = nodes_to_label_info[parent_cand * (cnt_unique_label + 1) + label_cur];
//                pair<int, int> query_result = nodesToLabelInfo(parent_cand, label_cur);
				for (int z = query_result.first; z < query_result.first + query_result.second; z++) {
					int can_id = nodes_data[z]; //v
					if (flag_prelin_char[can_id] == check_value){
						flag_prelin_char[can_id] ++;
						if (check_value == 0)
							array_to_clean[ to_clean_index ++ ] = can_id;
					}
				}
			}
			check_value ++; //update the check value by one
		}
		//end 1.1.

		//1.2. for each vertex v in u.C	
		NodeIndexUnit & cur_node_unit = indexSet[current_node];

		for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++){ 
			int can_id = cur_node_unit.candidates[cand_index]; //v
			if( flag_prelin_char[can_id] != check_value ) {//if v.cnt != Cnt then Remove v from u.C
				cur_node_unit.candidates[cand_index] = cur_node_unit.candidates[ --cur_node_unit.size];
				--cand_index;
			}
		}
		//end 1.2.

		//1.3. reset v.cnt = 0 for all v in G such that v.cnt > 0
		while(to_clean_index != 0) //reset v.cnt = 0 for all v in G s.t. v.cnt > 0;
			flag_prelin_char [ array_to_clean[ --to_clean_index ] ] = 0;
		//end 1.3.

	}//for simulation sequence
}

inline void sortByDeg(){
	for (int u = 0; u < cnt_node_query; ++u){
		NodeIndexUnit& cur_node_unit = indexSet[u];
		sort(cur_node_unit.candidates, cur_node_unit.candidates + cur_node_unit.size,
		[](int x, int y) -> bool
		{
#ifdef DEGORDER
			return degree_data[x] < degree_data[y];
#else
			return core_number_data[x] < core_number_data[y];
#endif
		});
	}
}
	
inline void adjacencyListConstruction (){
	//0. for each query vertex u initialize NodeIndexUnit
	for(int i = 0; i < simulation_sequence_index; i++)
	{
		int current_node = simulation_sequence[i]; // u
		int label_cur = nodes_label_query[current_node];

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node) continue;
		NodeIndexUnit& cur_node_unit = indexSet[current_node];

		cur_node_unit.backtrack_neighbor_cand_size = new int[node_degree_query[current_node]]();
		cur_node_unit.backtrack_size_of_index = new int*[node_degree_query[current_node]];
		cur_node_unit.backtrack_index = new int**[node_degree_query[current_node]];
		for(int index = 0; index < node_degree_query[current_node]; index++)
		{
			cur_node_unit.backtrack_size_of_index[index] = NULL;
			cur_node_unit.backtrack_index[index] = NULL;
		}
	}

	//1. for each query vertex u
	for(int i = 0; i < simulation_sequence_index; i++)
	{
		int current_node = simulation_sequence[i]; // u
		int label_cur = nodes_label_query[current_node];

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node) continue;
		NodeIndexUnit& cur_node_unit = indexSet[current_node];

		// cur_node_unit.backtrack_neighbor_cand_size = new int[node_degree_query[current_node]]();
		// cur_node_unit.backtrack_size_of_index = new int*[node_degree_query[current_node]];
		// cur_node_unit.backtrack_index = new int**[node_degree_query[current_node]];
		// for(int index = 0; index < node_degree_query[current_node]; index++)
		// {
		// 	cur_node_unit.backtrack_size_of_index[index] = NULL;
		// 	cur_node_unit.backtrack_index[index] = NULL;
		// }

		//1.0. preprocess u.C 1) to check if v in u.C in constant time 2) to know v's cand_index
		for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++)
		{
			int can_id = cur_node_unit.candidates[cand_index]; //v
			flag_child_cand[can_id] = cand_index;
		}

		//1.1. for each neighbor un of u
		for(int neighbor_index = 0; neighbor_index < node_degree_query[current_node]; neighbor_index++)
		{
			int neighbor = query_nodes_array[query_nodes_array_info[current_node] + neighbor_index];
			NodeIndexUnit& neighbor_node_unit = indexSet[neighbor];
			int u_index = query_neighbor_index[neighbor][current_node];
			
			// NEC boost
			if(NEC_map[neighbor] != -1 && NEC_map[neighbor] != neighbor) continue;

			memset(index_array_for_indexSet, 0, sizeof(int) * cur_node_unit.size);

			//1.1.0. initialize bactrack index
			if(neighbor_node_unit.backtrack_neighbor_cand_size[u_index] < cur_node_unit.size)
			{
				neighbor_node_unit.backtrack_size_of_index[u_index] = new int[cur_node_unit.size]();
				neighbor_node_unit.backtrack_index[u_index] = new int*[cur_node_unit.size];
				for(int index = 0; index < cur_node_unit.size; index++)
					neighbor_node_unit.backtrack_index[u_index][index] = NULL;

				neighbor_node_unit.backtrack_neighbor_cand_size[u_index] = cur_node_unit.size;
			}

			//compute index size for each v in u.C
			for(int neighbor_cand_index = 0; neighbor_cand_index < neighbor_node_unit.size; neighbor_cand_index++)
			{
				int neighbor_cand = neighbor_node_unit.candidates[neighbor_cand_index];
				pair<int, int> query_result = nodes_to_label_info[neighbor_cand * (cnt_unique_label + 1) + label_cur]; // information about vertex such that it is a neighbor of 'neighbor_cand' and its label is label_cur
				for(int i = query_result.first; i < query_result.first + query_result.second; i++)
				{
					int can_id = nodes_data[i];
					if(flag_child_cand[can_id] != -1) //if v in u.C
					{
						index_array_for_indexSet[flag_child_cand[can_id]]++;
					}
				}
			}

			//initialize bactrack_index for each v in u.C
			for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++)
			{
				if(neighbor_node_unit.backtrack_size_of_index[u_index][cand_index] < index_array_for_indexSet[cand_index])
				{
					neighbor_node_unit.backtrack_index[u_index][cand_index] = new int[index_array_for_indexSet[cand_index]];
				}
				neighbor_node_unit.backtrack_size_of_index[u_index][cand_index] = index_array_for_indexSet[cand_index];
			}
			memset(index_array_for_indexSet, 0, sizeof(int) * cur_node_unit.size);
			//end 1.1.0.

			//1.1.1. for each vn in un.C
			for(int neighbor_cand_index = 0; neighbor_cand_index < neighbor_node_unit.size; neighbor_cand_index++)
			{
				int neighbor_cand = neighbor_node_unit.candidates[neighbor_cand_index];
				
				//1.1.1.1. for each v adjacent to vn in G with label
				pair<int, int> query_result = nodes_to_label_info[neighbor_cand * (cnt_unique_label + 1) + label_cur];
				//for each vertex v in N_G(vn) with label l_q(u)
				for(int z = query_result.first; z < query_result.first + query_result.second; z++)
				{
					int can_id = nodes_data[z]; //v
					if(flag_child_cand[can_id] != -1) //if v in u.C
					{
						int cand_index = flag_child_cand[can_id];
						neighbor_node_unit.backtrack_index[u_index][cand_index][index_array_for_indexSet[cand_index]++] = neighbor_cand_index;
					}
				}
			}

			for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++)
			{
				neighbor_node_unit.backtrack_size_of_index[u_index][cand_index] = index_array_for_indexSet[cand_index];
			}
		}

		//1.2. reset v in G such that v in u.C and v.index = v's index in u.C
		for(int cand_index = 0; cand_index < cur_node_unit.size; cand_index++)
		{
			int can_id = cur_node_unit.candidates[cand_index];
			flag_child_cand[can_id] = -1;
		}

	}
}

void clearMemory()
{
	//cout << "clearing memory... root node id=" << root_node_id << endl<< flush;
	for(int i = 0; i < cnt_node_query; i++)
	{
		NodeIndexUnit& current_node = indexSet[i];
		//cout << "clear " << i <<"th node " << endl << flush; 
		if(current_node.backtrack_index != NULL)
		{
			for(int neighbor_index = 0; neighbor_index < node_degree_query[i]; neighbor_index++)
			{
				if(current_node.backtrack_index[neighbor_index] != NULL)
				{
					for(int cand_index = 0; cand_index < current_node.backtrack_neighbor_cand_size[neighbor_index]; cand_index++)
					{
						if(current_node.backtrack_index[neighbor_index][cand_index] != NULL)
							delete[] current_node.backtrack_index[neighbor_index][cand_index];
							current_node.backtrack_index[neighbor_index][cand_index] = NULL;
					}
					delete[] current_node.backtrack_size_of_index[neighbor_index];
					current_node.backtrack_size_of_index[neighbor_index] = NULL;
					delete[] current_node.backtrack_index[neighbor_index];
					current_node.backtrack_index[neighbor_index] = NULL;
				}
			}

			if(current_node.backtrack_index != NULL)
			{
				delete[] current_node.backtrack_neighbor_cand_size;
				current_node.backtrack_neighbor_cand_size = NULL;
				delete[] current_node.backtrack_index;
				current_node.backtrack_index = NULL;
			}

		}
	}
    memset(chosen_vertex, false, sizeof(bool) * cnt_node);
    /*
    if(nodes_label != NULL) { delete[] nodes_label; nodes_label = NULL; }
    if(label_freqency != NULL) { delete[] label_freqency; label_freqency = NULL; }
    if(transferredLabel != NULL) { delete[] transferredLabel; transferredLabel = NULL; }
    if(nodes_info != NULL) { delete[] nodes_info; nodes_info = NULL; }
    if(nodes_to_label_info != NULL) { delete[] nodes_to_label_info; nodes_to_label_info = NULL; }
    if(core_number_data != NULL) { delete[] core_number_data; core_number_data = NULL; }
    if(degree_data != NULL) { delete[] degree_data; degree_data = NULL; }
    if(nodes_data != NULL) { delete[] nodes_data; nodes_data = NULL; }
    if(NLF_check != NULL) { delete[] NLF_check; NLF_check = NULL; }
    if(NLF_array != NULL) { delete[] NLF_array; NLF_array = NULL; }
    if(MAX_NB_degree_data != NULL) { delete[] MAX_NB_degree_data; MAX_NB_degree_data = NULL; }
    if(NEC_mapping != NULL) { delete[] NEC_mapping; NEC_mapping = NULL; }
    if(NEC_mapping_Actual != NULL) { delete[] NEC_mapping_Actual; NEC_mapping_Actual = NULL; }
    if(NEC_mapping_pair != NULL) { delete[] NEC_mapping_pair; NEC_mapping_pair = NULL; }
    if(NEC_set_array != NULL) { delete[] NEC_set_array; NEC_set_array = NULL; }
    */
	//cout << "done!" << endl << flush;
}

void buildGlobalMatchingOrder()
{
	memset(visited_for_query, 0, sizeof(char) * MAX_QUERY_NODE);
	int min_C_index = -1;
	for(int i = 0; i < cnt_node_query; i++)
	{
		if(NEC_map[i] != -1) continue;
		if(min_C_index == -1) min_C_index = i;
		else if(indexSet[i].size < indexSet[min_C_index].size) min_C_index = i;
	}
	
	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
	pq.push(make_pair(indexSet[min_C_index].size, min_C_index));
	visited_for_query[min_C_index] = true;
	int idx = 0;
	while(!pq.empty())
	{
		int u = pq.top().second;
		pq.pop();
		for(int i = query_nodes_array_info[u]; i < query_nodes_array_info[u + 1]; i++)
		{
			int v = query_nodes_array[i];
			if(!visited_for_query[v] && NEC_map[v] == -1)
			{
				visited_for_query[v] = true;
				pq.push(make_pair(indexSet[v].size, v));
			}
		}
		global_matching_order[idx++] = u;
	}
}

//==========================================================================================

inline bool BFS_MM(){ //BFS function for maximum matching

	/*
	 * Note that in this function, we are NOT manipulating the actual query nodes and data nodes.
	 * Instead, we are dealing with the index in LLC and the data node index in the sum NEC candidate set.
	 */

	int * queue = global_temp_array_1;
	int queue_start = 0;//use for popping up
	int queue_end = 0;//use for pushing back

	for (int u = 1; u < nodes_MM_idx; u++){
		if ( pair_U[u] == NIL ){
			dist[u] = 0;
			queue[queue_end ++] = u;
		}
		else
			dist[u] = INF;
	}

	dist[NIL] = INF;

	while (queue_start != queue_end){

		int u = queue[queue_start];
		queue_start ++;

		if (dist[u] < dist[NIL]){

			int nec = nec_mappiing_MM[u]; //get the actual index in the LLC

			for (int j = 0; j < nec_region_size[nec]; j++){
				int v =  u_cand_MM[nec * max_label_counter + j];
				if (dist [ pair_V[v] ] == INF){
					dist[ pair_V[v] ] = dist[ u ] + 1; //here use "u" instead of "nec", because "nec" is only used to get the candidate set.
					queue[queue_end ++] = pair_V[v];
				}
			}
		}
	}

	return (dist[ NIL ] != INF);
}

inline bool DFS_MM(int u){

	if (u != NIL){

		int nec = nec_mappiing_MM[u]; //get the actual index in the LLC

		for (int j = 0; j < nec_region_size[nec]; j++){
			int v =  u_cand_MM[nec * max_label_counter + j];
			if (dist[ pair_V[ v ] ] == dist[ u ] + 1)
				if (DFS_MM( pair_V[ v ] )){
					pair_V[ v ] = u;
					pair_U[ u ] = v;
					return true;
				}
		}

		dist[ u ] = INF;
		return false;
	}

	return true;
}

inline void combine( double & count_mapping, double result_so_far_from_former_llc){

	/*
	 * max_cand version of combination
	 */
	
	if (LIMIT > 0 and count_mapping * result_so_far_from_former_llc > LIMIT_REMAINED)
		return;

	//Step One: find a candidate that is with the max number of NEC nodes' candidate set containing it
	//locate the candidate with the last coverage number
	int max_cand = -1;
	int max_cand_size = 1;
	int max_cand_idx = -1;
	int * local_sum_nec_cand = sum_nec_cands_array[curretn_LLC_label];
	int & local_size_sum = size_sum_nec_cands[curretn_LLC_label];

	for (int i = 1; i < local_size_sum; i++){
		int can = local_sum_nec_cand[i];
		int temp = idx_v_cands[ can ];
		if ((int) mapping_flag_data[can] == 0 && temp > max_cand_size){
			max_cand_size = temp;
			max_cand_idx = i;
		}
	}

	if (max_cand_idx == -1){
		//in this case, there is no such a cand that has more than one query node having it as a candidate
		//hence, we can perform the combination and permutation to get the result here
		double res = 1;
		for (int i = 0; i < nec_count_set_size; i++){
			if (nec_count_set[i] != local_flag_query[i]){//this node is unmatched or still has unmatched NEC nodes
				int rest = nec_count_set[i] - local_flag_query[i];
				int remaining_nec_region_size = nec_region_size[i];
				if (remaining_nec_region_size < rest) //smaller ==> no result
					return;
				if (remaining_nec_region_size == rest) //equal ==> only one result
					res *= 1;
				else {	//larger than rest, then perform the combination
					for (int x = remaining_nec_region_size - rest + 1; x <= remaining_nec_region_size; x++)
						res *= x;
					res /= factorization(rest);
				} //end else
			}
		}

		res *= PRE_COMPUTED_PERMUTATION;
		count_mapping += res;
		return;
	} //end if max_cand = -1

	max_cand = local_sum_nec_cand[max_cand_idx];//the case that max_cand is not -1
	swap_value(local_sum_nec_cand[local_size_sum - 1], local_sum_nec_cand[max_cand_idx]);//swap to the end
	local_size_sum --;	//reduce the size by one

	//then we examine the nec nodes that are covered by the max cand
	/*
	 * The pruning rule here is that :
	 * 	for this cand, if there exists multiple query nodes covered by it and those nodes' cand size (the conditional cand size for the nec node) is 1,
	 * 	then there is no result. Because, one cand cannot be mapped to multiple query node at the same time.
	 * 	But if there only exists one such node, then it is OK to precede and this cand should be mapped to this query node.
	 */

	//now we try to map this cand node, and remove it from all candidate regions that contain it

	//this for loop remove the max_cand from the candidate sets that contain it
	for (int i = 0; i < idx_v_cands[max_cand]; i++){//for each query node whose cand set containing it
		//q_node_index is the position of this max_cand in the nec node set
		int q_node_index = v_cands_info[max_cand][i].first;
		int pos_in_region = v_cands_info[max_cand][i].second;
		int size_of_region = nec_region_size[q_node_index];
		if (pos_in_region ==  size_of_region - 1)
			nec_region_size[q_node_index] --;
		else if (pos_in_region <  size_of_region - 1){
			pair<int, int> & temp_p = u_cand_info [q_node_index * max_label_counter +  size_of_region - 1 ];
			int last_cand = temp_p.first;
			int last_pos = temp_p.second;
			swap_value( temp_p, u_cand_info [q_node_index * max_label_counter + pos_in_region ] );
			v_cands_info[max_cand][i].second = size_of_region - 1;
			v_cands_info[last_cand][last_pos].second = pos_in_region ;
			nec_region_size[q_node_index] --;
		}
	} // end for i

	//======multiple=====
	//sort accord to size of each query node's effective candidate size
	sort (v_cands_info[max_cand].begin(), v_cands_info[max_cand].begin() + idx_v_cands[max_cand], sort_by_effective_cand_size_pair);

	//the following loop maintains the position info after the sorting
	for (int i = 0; i < idx_v_cands[max_cand]; i++){ //for each query node containing max_cand
		pair<int, int> temp_p =  v_cands_info[max_cand][i];
		int q_node_index = temp_p.first;
		int pos_in_region =  temp_p.second;
		u_cand_info[q_node_index * max_label_counter + pos_in_region].second = i;//set the position of "q_node_index" in the u_cand set of the candidate "max_cand"
	}// end for i


	for (int i = 0; i < idx_v_cands[max_cand]; i++){ //for each query node containing max_cand
		int q_node_index = v_cands_info[max_cand][i].first;
		if ((int)local_flag_query[q_node_index] == nec_count_set[q_node_index])//check whether a nec node has been fully mapped or not
			continue;
		mapping_flag_data[max_cand] = 1; //map "max_cand" to this query node
		local_flag_query[q_node_index] ++; //increase the flag value accordingly

		//here, if the query node "q_node_index" has been fully mapped, then we need to delete it from all the data node's candidate set which
		//containing it, such  that this query node will not be counted in the future round of finding the next max_cand
		if (local_flag_query[q_node_index] == nec_count_set[q_node_index]){//this query node is fully mapped

			for (int j = 0; j < nec_region_size[q_node_index]; j++){//for each of its candidate
				int v_node = u_cand_info[q_node_index * max_label_counter + j].first; // get the candidate from position "j"
				vector < pair<int, int>> & v_c_info = v_cands_info[v_node];
				int pos = u_cand_info[q_node_index * max_label_counter + j].second;
				if (pos == idx_v_cands[v_node] - 1)
					idx_v_cands[ v_node ] --;
				else if (pos < idx_v_cands[v_node] - 1){
					int size = idx_v_cands[v_node];
					swap_value( u_cand_info[ v_c_info [ size - 1].first * max_label_counter +  v_c_info [ size - 1 ].second ],
							u_cand_info[ v_c_info [ pos ].first * max_label_counter + v_c_info [ pos ].second] );	//auxiliary operation
					swap_value ( v_c_info [ size - 1 ], v_c_info [ pos ]);	//auxiliary operation
					idx_v_cands[ v_node ] --;
				}
			}// end for j
		}// end if

		combine(count_mapping, result_so_far_from_former_llc); //recursive function

		if (LIMIT > 0 and count_mapping * result_so_far_from_former_llc > LIMIT_REMAINED) {//clean up and return
			mapping_flag_data[max_cand] = 0;
			return;
		}

		//now we have to update the c(v) again, by putting back the last mapped query node (actually here: the q_node_index)
		if (local_flag_query[q_node_index] == nec_count_set[q_node_index]){//NEC equal
			for (int j = 0; j < nec_region_size[q_node_index]; j++){
				int v_node = u_cand_info[q_node_index * max_label_counter + j].first;
				if (v_node == max_cand) //skip "max_cand"
					continue;
				idx_v_cands[ v_node ] ++; //extend the size by one
			}
		}
		local_flag_query[q_node_index] --;//reduce this query node's flag, because it has been unmapped.
	}

	mapping_flag_data[max_cand] = 0;
	combine(count_mapping, result_so_far_from_former_llc); //recursive function

	if (LIMIT > 0 and count_mapping * result_so_far_from_former_llc > LIMIT_REMAINED) {//clean up and return
		mapping_flag_data[max_cand] = 0;
		return;
	}
	// end =============== (multiple)

	//put back the max cand into each query node's candidate set
	for (int i = 0; i < idx_v_cands[max_cand]; i++){
		int q_node_index = v_cands_info[max_cand][i].first;
		nec_region_size[q_node_index] ++;
	}

	local_size_sum ++;
	mapping_flag_data[max_cand] = 0;
}

inline void init_weight_by_min()
{
	pair<int, int> pos;
  	int start, end, v, vc_index, u, uc;
	long long min, total;
	bool flag = false;
  	for(int level = level_index.size() - 1; level >= 0; level--){
    	pos = level_index[level];
	    start = pos.first;
	    end = pos.second;
	    //for(int seq_index = start; seq_index < end; seq_index++){
	    for(int seq_index = end - 1; seq_index >= start; seq_index--){
    		u = simulation_sequence[seq_index];
		    if(NEC_map[u] != -1)
        		continue;
	        NodeIndexUnit & cur_node_unit = indexSet[u];
      		for(int i = 0; i < cur_node_unit.size; ++i){ // for each vertex v in u.C
		        v = cur_node_unit.candidates[i];     // v = u.C
		        if(v < 0) continue;
		        //cout<<"v: "<<v<<endl;
		        min = LLONG_MAX;
		        flag = false;
		        for(int j = 0; j < DAG_child_query_size[u]; ++j){ // for each child uc of u in dag(q)
			        uc = DAG_child_query[u][j];
			        if(NEC_map[uc] != -1)
				        continue;
					//cout<<"uc: "<<uc<<", DAG_parent_query_size: "<<DAG_parent_query_size[uc]<<endl;
			        if(DAG_parent_query_size[uc] > 1 ) continue;
		        	flag = true;
		            total = 0;
			        NodeIndexUnit & child_node_unit = indexSet[uc];
					int u_index = query_neighbor_index[uc][u];
	        	    //cout << "|N^{"<<u<<"}_{"<<uc<<"}("<<v<<")|: "<<child_node_unit.backtrack_size_of_index[0][i]<<endl;
		            for(int k = 0; k < child_node_unit.backtrack_size_of_index[u_index][i]; ++k){ // for each vc in N^{u}_{uc}(v)
			            vc_index = child_node_unit.backtrack_index[u_index][i][k];
			            if(vc_index < 0) continue;
		            	//cout << "N^{"<<u<<"}_{"<<uc<<"}("<<v<<")["<<k<<"]: "<<child_node_unit.candidates[vc_index]<<endl;
			            if (child_node_unit.weight[vc_index] > LLONG_MAX-1 - total){
	            			//cout<<"error"<<endl;
					        //exit(1);

							total = LLONG_MAX-1;
							break;
		    			}
			            total += child_node_unit.weight[vc_index];
	    	        	//cout << "W_{"<<uc<<"}("<<vc_index<<"): "<<child_node_unit.weight[vc_index]<<endl;
        			}
		        	if(total < min) min = total;
				}
		        if(!flag) 
			        cur_node_unit.weight[i] = 1;
		        else 
        			cur_node_unit.weight[i] = min;
		        //cout << "W_"<<u<<"("<<v<<") = "<<cur_node_unit.weight[i]<<endl;
			}
		}
	}  
}

inline double matchLeaves(double already_found_mapping) {
	/*
	 *  with maximum matching and max_cand optimization.
	 */

	LIMIT_REMAINED = LIMIT - already_found_mapping;

	memset(nec_mappiing_MM, -1, sizeof (int) * (cnt_node_query + 1) ); //for MM
	double actual_mapping_result = 1;
	int nec_set_size = NEC_set_by_label_index.size() - 1;// -1 is because the last element is actually redundant
	int can; //a candidate
	int unique_u_cand_counter;
	int sum_nec_num; // the sum number of all nec node of a label
	int u_cand_counter;
	int sum_cand;

	Leaf_cands_idx = 0;

	//=== First scan === fastly identify whether all candidates satisfy the (cand < nec) and (sum_cand < sum_nec). (the sum cands here are not unique cands)
	for (int i = 0; i < nec_set_size; i++) { //for each label

		int start = NEC_set_by_label_index[i].second; // start position
		int end = NEC_set_by_label_index[i + 1].second; // end position
		sum_nec_num = 0;
		sum_cand = 0;
		for (int j = start; j < end; j++) {
			//basic info about this nec node
			int parent_id = NEC_set_array[j].parent_id;
			int nec_num = NEC_set_array[j].sum;
			int represent_node = NEC_set_array[j].represent_node;
			sum_nec_num += nec_num;
			u_cand_counter = 0; //record the number of candidate for this nec node
			int parent_pos = self_pos[parent_id];
			NodeIndexUnit & unit = indexSet[represent_node];
			for (int it = 0; it < unit.backtrack_size_of_index[0][parent_pos]; it++){ // 20170414
                int v_id = unit.candidates[unit.backtrack_index[0][parent_pos][it]]; // 20170414
                //XXX
				if (v_id >= 0 and mapping_flag_data[ v_id ] == 0) // 20170414
					u_cand_counter++;
            }
			if (u_cand_counter < nec_num)
				return 0;
			sum_cand += u_cand_counter;
		}
		if (sum_cand < sum_nec_num)
			return 0;
	}

	//=== Second scan: extract cands and perform the maximum matching
	for (int i = 0; i < nec_set_size; i++) { //for each label

		int start = NEC_set_by_label_index[i].second; // start position
		int end = NEC_set_by_label_index[i + 1].second; // end position
		int label = NEC_set_by_label_index[i].first; //the label of this nec set

		//initialization
		idx_sum_nec_cands = 1;
		//flag_sum_nec_cands must be all -1 here.
		unique_u_cand_counter = 0;
		sum_nec_num = 0;
		sum_cand = 0;
		nodes_MM_idx = 1; //for MM

		int * local_sum_nec_cand = sum_nec_cands_array[label];

		for (int j = start; j < end; j++) {

			//basic info about this nec node
			int parent_id = NEC_set_array[j].parent_id;
			int nec_num = NEC_set_array[j].sum;
			int represent_node = NEC_set_array[j].represent_node;
			sum_nec_num += nec_num;
			for (int mm = 0; mm < nec_num; mm++) //for MM
				nec_mappiing_MM [nodes_MM_idx ++] = j;
			u_cand_counter = 0; //record the number of candidate for this nec node
			int parent_pos = self_pos[parent_id];
			NodeIndexUnit & unit = indexSet[represent_node];
			int cand_start = Leaf_cands_idx;

			for (int it = 0; it < unit.backtrack_size_of_index[0][parent_pos]; it++){ // 20170414
				int can = unit.candidates[ unit.backtrack_index[0][parent_pos][it] ]; // 20170414
				if (mapping_flag_data[ can ] == 0){
					if (flag_sum_nec_cands[can] == -1){ //if not already in the candidate set,
						flag_sum_nec_cands[can] = idx_sum_nec_cands; //and mark it true and indicate its position in "sum_nec_cands"
						local_sum_nec_cand[idx_sum_nec_cands ++] = can; //then add it into the candidate's "idx_sum_nec_cands" position,
						unique_u_cand_counter ++;	//optimization added on 21:10 05/05/2016
					}
					u_cand_MM[j * max_label_counter + u_cand_counter] = flag_sum_nec_cands[can]; //for the maximum matching
					u_cand_counter ++;
					Leaf_cands[Leaf_cands_idx ++] = can; //store the candidate for the second stage (combination)
				}
			}

			Leaf_cands_info[ j ] = make_pair(cand_start, u_cand_counter);
			nec_region_size[j] = u_cand_counter;	//set the size of the j-th candidate set of this nec node
			sum_cand += u_cand_counter;
			size_sum_nec_cands[label] = idx_sum_nec_cands; //added on 20:15 09/05/2016
		}

		if (unique_u_cand_counter < sum_nec_num){//no result //reset "flag_sum_nec_cands"
			for(int i = 1; i < idx_sum_nec_cands; i++) //it starts from 1 NOT 0
				flag_sum_nec_cands[ local_sum_nec_cand[i] ] = -1;
			return 0;
		}

		//after the  "sum_nec_cands" is initialized, we reset the "flag_sum_nec_cands" for next time use
		for(int i = 1; i < idx_sum_nec_cands; i++) //it starts from 1 NOT 0
			flag_sum_nec_cands[ local_sum_nec_cand[i] ] = -1;

		//== BEGIN checking MAXIMUM MATCHING ==============================================
		//using nec_region as the adjacent list
		memset(pair_U, NIL, sizeof (int) * ( cnt_node_query + 1) ); //must reset before using
		memset(pair_V, NIL, sizeof (int) * max_label_counter); //must reset before using

		int matching = 0;

		while ( BFS_MM() ){
			for (int u = 1; u < nodes_MM_idx; u++)
				if (pair_U[u] == NIL)
					if ( DFS_MM(u) )
						matching ++;
		}

		if (matching != nodes_MM_idx - 1)
			return 0;
		//== END checking MAXIMUM MATCHING ==========================================
		NEC_set_ranking[i] = make_pair(i, make_pair(sum_nec_num, sum_cand));
	}

	sort(NEC_set_ranking, NEC_set_ranking + nec_set_size, sortLLCnode);

	//Deal with the labels one by one
	for (int i = 0; i < nec_set_size; i++) { //for each nec node label

		if (LIMIT > 0 and actual_mapping_result > LIMIT_REMAINED)//dont continue to find matchings, because we know there would be at least one matching for this label due to MM
			continue;

		double count_local_mapping_label = 0;
		int label_index = NEC_set_ranking[i].first;	// the label index of this nec label set ==> this is not the actual label!!!!
		int node_sum = NEC_set_ranking[i].second.first; // the number of node in this nec label set
		curretn_LLC_label = NEC_set_by_label_index[label_index].first;

		if (node_sum == 1) {	//==== CASE ONE : there is only one node in this nec set with this label =====
			count_local_mapping_label = (long long)NEC_set_ranking[i].second.second; //we have computed this one in the last step
			if (count_local_mapping_label == 0)
				return 0;
		} else {	//==== CASE TWO : more than one node, and possible more than one start nodes (nec_units)

			memset(local_flag_query, 0, sizeof(char) * cnt_node_query);
			idx_sum_nec_cands = 1;//omit zero, start from 1 ===> Leave 0 for NIL of Maximum Matching
			int start = NEC_set_by_label_index[label_index].second;
			int end = NEC_set_by_label_index[label_index + 1].second;
			int nec_size = end - start; //number of nec this label has

			//======= sorting here is to decide the processing sequence of the nec nodes in this nec label set
			for (int j = start, x = 0; j < end; j++, x++)
				v_nec_count[x] = make_pair(j, NEC_set_array[j].sum);
			sort(v_nec_count, v_nec_count + nec_size, sort_by_second_element);

			int * to_clean_idx_v_cands = global_temp_array_2; //store the candidates that need to be cleaned for the array "idx_v_cands"
			int idx_to_clean_idx_v_cands = 0;
			//in this loop, for this point beyond, before any return, "idx_v_cands" must be cleaned up using the above two parameters

			nec_count_set_size = nec_size;

			for (int j = 0; j < nec_size; j++) {

				int nec_index = v_nec_count[j].first;
				int nec_count = v_nec_count[j].second;
				nec_count_set[j] = nec_count; // the sum of nodes that the nec representative stands for
				u_cand_counter = 0;

				pair <int, int> temp_p = Leaf_cands_info[ nec_index ];
				for(int x = temp_p.first; x < temp_p.first + temp_p.second; x++){

					can = Leaf_cands[x];
					int temp = idx_v_cands[can]; //get the position for storing this candidate

					//the temp-th position of v_cands_info[can] stores a candidate "can" which is for "j-th" nec node and
					//stored in the "u_cand_counter-th" position in the candidate set
					v_cands_info[can][temp] = make_pair(j, u_cand_counter); // pair<j-th nec node, can's pos in the candidate set of j-th nec node>

					//store the candidate and its corresponding info ==>  <the candidate can, the candidate's pos in idx_v_cands[can]>
					//in the single-array: "j-th" nec node and the "u_cand_counter-th" position.
					u_cand_info[j * max_label_counter + u_cand_counter] = make_pair(can, temp);

					if(temp == 0)//record for the cleaning purpose
						to_clean_idx_v_cands[ idx_to_clean_idx_v_cands ++] = can;

					idx_v_cands[can] ++; // update the counter(next available position) for "idx_v_cands[can]"
					u_cand_counter ++; //update the candidate counter(next available position)
				}

				//if the number of candidate for this nec node, is smaller than the nec count of this nec node, then no result can be found.
				nec_region_size[j] = u_cand_counter;
			}//end for j


			//to reduce the computation cost, we pre-compute the value of the all factorization
			PRE_COMPUTED_PERMUTATION = 1;
			for (int i = 0; i < nec_count_set_size; i++)
				PRE_COMPUTED_PERMUTATION *= factorization(nec_count_set[i]);

			combination_array_index = 0;

			int * to_clean_mapping_flag_data = global_temp_array_1; //store the candidates that have been set to mapped in "mapping_flag_data"; and to be cleaned before returning
			int idx_to_clean_mapping_flag_data = 0;
			//in this loop, for this point beyond, before any return, "mapping_flag_data" must be cleaned up using the above two parameters.

			int found = 1;
			while (found > 0){
				//=============preprocess before the combination
				found = 0;

				//first, we scan each nec node, to see whether there exists such a node that its candidate size is equal to its unmatched nec size
				for (int i = 0; i < nec_count_set_size; i++){

					if (nec_region_size[i] != 0 && nec_count_set[i] > nec_region_size[i]){
						while (idx_to_clean_mapping_flag_data != 0) //cleaning up mapping_flag_data that changed during preprocessing
							mapping_flag_data[ global_temp_array_1[ -- idx_to_clean_mapping_flag_data ] ] = 0;
						while (idx_to_clean_idx_v_cands != 0) //cleaning up idx_v_cands
							idx_v_cands[ to_clean_idx_v_cands[ -- idx_to_clean_idx_v_cands ] ] = 0;
						return 0;
					}

					//all can be matched one on one
					if (nec_count_set[i] == nec_region_size[i]){
						found ++;
						local_flag_query[i] = nec_count_set[i];//mark this nec node all matched

						//update the c(v); there is no neccessary to update the c(v)
						for (int j = 0; j < nec_region_size[i]; j++){

							int v_node = u_cand_info[i * max_label_counter + j].first; //the j-th candidate
							mapping_flag_data[v_node] = 1;//set this cand's flag to 1
							to_clean_mapping_flag_data[idx_to_clean_mapping_flag_data ++] = v_node;//recording this cand

							//for each query node that contains this candidate "v_node"
							for (int x = 0; x < idx_v_cands[v_node]; x++){ // x is the position in "v_cands_info[v_node]"

								int q_node_index = v_cands_info[v_node][x].first;//get the original index in the nec node set
								if (q_node_index == i)
									continue;

								//the position of v_node in the corresponding("q_node_index-th" ) candidate set.
								int pos_in_region = v_cands_info[v_node][x].second;
								int size_of_region = nec_region_size[q_node_index];

								if (pos_in_region ==  size_of_region - 1) //this is the last node in the candidate set
									nec_region_size[q_node_index] --; //remove this candidate by downsizing the nec_region_size by one
								else if (pos_in_region <  size_of_region - 1){
									//normal case: not the last one.
									//swap the two candidates and their corresponding info, then downsize the nec_region_size by one
									pair<int, int> & temp_p = u_cand_info [q_node_index * max_label_counter +  size_of_region - 1 ];
									int last_cand = temp_p.first;
									int last_pos = temp_p.second;
									swap_value ( temp_p, u_cand_info [ q_node_index * max_label_counter + pos_in_region ]); //swap the candidate
									v_cands_info[v_node][x].second = size_of_region - 1;
									v_cands_info[last_cand][last_pos].second = pos_in_region;
									nec_region_size[q_node_index] --;
								}


								if (nec_region_size[q_node_index] < nec_count_set[q_node_index]){
									while (idx_to_clean_mapping_flag_data != 0) //cleaning up mapping_flag_data that changed during preprocessing
										mapping_flag_data[ to_clean_mapping_flag_data[ -- idx_to_clean_mapping_flag_data ] ] = 0;
									while (idx_to_clean_idx_v_cands != 0) //cleaning up idx_v_cands
										idx_v_cands[ to_clean_idx_v_cands[ -- idx_to_clean_idx_v_cands ] ] = 0;
									return 0;
								}
							} // end for x
						}// end for j
						nec_region_size[i] = 0;
					}// end if
				}// end for i
			}//end PREPROCESSING
			//================== END PREPROCESSING ==============================

			combine(count_local_mapping_label, actual_mapping_result);

			//================= cleaning up =========================================
			while (idx_to_clean_mapping_flag_data != 0) //cleaning up mapping_flag_data that changed during preprocessing
				mapping_flag_data[ to_clean_mapping_flag_data[ -- idx_to_clean_mapping_flag_data ] ] = 0;
			while (idx_to_clean_idx_v_cands != 0) //cleaning up idx_v_cands
				idx_v_cands[ to_clean_idx_v_cands[ -- idx_to_clean_idx_v_cands ] ] = 0;
			//=============================================================================

			if (count_local_mapping_label == 0)
				return 0;

		} //end else : end of case TWO

		actual_mapping_result *= count_local_mapping_label;

		if (LIMIT > 0 and actual_mapping_result >= LIMIT_REMAINED)
			continue;
	}
	return actual_mapping_result;
}

inline void reinsert_to_frontier_improved(int u, weight_type weight, int position)
{
	if(frontier_node_idx != position)
	{
		frontier_node[frontier_node_idx] = frontier_node[position];
		weight_array[frontier_node_idx] = weight_array[position];
		frontier_index[frontier_node[frontier_node_idx]] = frontier_node_idx;
	}
    frontier_node_idx++;

    frontier_node[position] = u;
    weight_array[position] = weight;
	frontier_index[u] = position;
    if(use_path_size_order)
    {
        frontier_min_index = position;
        frontier_min_weight = weight;
    }
    #ifdef MAPPING_FUNCTION_LOG
    cout<<"(I)heapSize: "<<frontier_node_idx<<endl;
    assert(frontier_node_idx >= 0);
    #endif
}

inline void insert_to_frontier_first_improved(int u, weight_type weight)
{
    if(use_path_size_order)
    {
        if (DAG_parent_query_size[u] > 1){
            frontier_min_index = frontier_min_weight = -1;
        }
        else if (weight < frontier_min_weight) {
            frontier_min_weight = weight;
            frontier_min_index = frontier_node_idx;
        }
    }
    frontier_node[frontier_node_idx] = u;
    weight_array[frontier_node_idx] = weight;
	frontier_index[u] = frontier_node_idx;
    frontier_node_idx++;
    #ifdef MAPPING_FUNCTION_LOG
    cout<<"(I)heapSize: "<<frontier_node_idx<<endl;
    assert(frontier_node_idx >= 0);
    #endif
}

inline void pop_from_frontier_improved(int& current, weight_type& cur_weight, int& position)
{
    int min_idx = -1;
    weight_type min_weight = WEIGHT_MAX;
    if(use_path_size_order)
    {
        if(frontier_min_index != -1){
            min_idx = frontier_min_index;
        } else {
            for(int i = 0; i < frontier_node_idx; i++){
                if(weight_array[i] < min_weight){
                    min_idx = i;
                    min_weight = weight_array[i];
                }
            }
        }
    }
    else if(use_candidate_size_order)
    {
        for(int i = 0; i < frontier_node_idx; i++){
            if(weight_array[i] < min_weight){
                min_idx = i;
                min_weight = weight_array[i];
            }
        }
    }

    current = frontier_node[min_idx];
    cur_weight = weight_array[min_idx];
    position = min_idx;
	frontier_index[frontier_node[min_idx]] = -1;

    frontier_node[min_idx] = frontier_node[frontier_node_idx - 1];
    weight_array[min_idx] = weight_array[frontier_node_idx - 1];
	frontier_index[frontier_node[min_idx]] = min_idx;
    frontier_node_idx--;
    #ifdef MAPPING_FUNCTION_LOG
    cout<<"(P)weight(u"<<current<<"): "<<cur_weight<<", heapSize: "<<frontier_node_idx<<endl;
    assert(frontier_node_idx >= 0);
    #endif
}

inline void remove_recently_inserted_from_frontier_improved(int n)
{
	// for(int i = 0; i < n; i++)
	// {
	// 	frontier_index[frontier_node[--frontier_node_idx]] = -1;
	// }
    frontier_node_idx -= n;
    #ifdef MAPPING_FUNCTION_LOG
    //cout << "C.7. Remove just inserted frontier elements "<< temp << endl;
    cout<<"(R)heapSize: "<<frontier_node_idx<<endl;
    assert(frontier_node_idx >= 0);
    #endif
}

inline void clear_frontier_improved()
{
    frontier_node_idx = 0;
    if(use_path_size_order)
    {
        frontier_min_index = -1;
        frontier_min_weight = LLONG_MAX;
    }
    #ifdef MAPPING_FUNCTION_LOG
    cout<<"(R)heapSize: "<<frontier_node_idx<<endl;
    assert(frontier_node_idx >= 0);
    #endif
}


inline long long compute_path_size_weight(int u, int num_mapped_parents, int up)
{
    long long weight = 0;
    int parent, parent_cand_index, parent2, parent_cand_index2;
    NodeIndexUnit & unit = indexSet[u];
    if(num_mapped_parents == 1){
		temp_intersection_array_index[u][0] = 0;
    	parent = DAG_parent_query[u][0];
    	parent_cand_index = self_pos[parent];
		int up_index = query_neighbor_index[u][parent];

        for(int j = 0; j < unit.backtrack_size_of_index[up_index][parent_cand_index]; ++j){ 
            int cand_idx = unit.backtrack_index[up_index][parent_cand_index][j];

			if(use_failing_set) 
				temp_intersection_array[u][0][temp_intersection_array_index[u][0]++] = cand_idx;
			if(!use_failing_set && !mapping_flag_data[unit.candidates[cand_idx]]) 
				temp_intersection_array[u][0][temp_intersection_array_index[u][0]++] = cand_idx;

            if(!mapping_flag_data[unit.candidates[cand_idx]]){
                if(weight > LLONG_MAX-1 -  unit.weight[cand_idx]){
                    //cout<<"overflow"<<endl;
                    //exit(1);
			        weight = LLONG_MAX-1;
			        break;
                }
                weight += unit.weight[cand_idx];
                //cout<<"[After] Case 1. weight += unit.weight["<<cand_idx<<"]("<<unit.weight[cand_idx]<<") -> weight: " << weight<<endl;
            }
        }
        //cout<<"[After] Case 1. return val: "<<weight<<endl;
        return weight;
    }
    else{
        temp_intersection_array_index[u][num_mapped_parents-1] = 0;
        parent = up;
        parent_cand_index = self_pos[parent];
        int up_index = query_neighbor_index[u][parent];

        int j = 0, k = 0, cand_idx, cand_idx2, backtrack_size, backtrack_size2;
        backtrack_size = temp_intersection_array_index[u][num_mapped_parents-2];
        backtrack_size2 = unit.backtrack_size_of_index[up_index][parent_cand_index];
        while(j < backtrack_size and k < backtrack_size2){
            cand_idx = temp_intersection_array[u][num_mapped_parents-2][j];
            cand_idx2 = unit.backtrack_index[up_index][parent_cand_index][k];
            if(cand_idx == cand_idx2){
				if(use_failing_set)
					temp_intersection_array[u][num_mapped_parents-1][temp_intersection_array_index[u][num_mapped_parents-1]++] = cand_idx;
				if(!use_failing_set && !mapping_flag_data[unit.candidates[cand_idx]])
					temp_intersection_array[u][num_mapped_parents-1][temp_intersection_array_index[u][num_mapped_parents-1]++] = cand_idx;

				if(!mapping_flag_data[unit.candidates[cand_idx]]){
					if(weight > LLONG_MAX-1 -  unit.weight[cand_idx]){
						//cout<<"overflow"<<endl;
						//exit(1);
						weight = LLONG_MAX-1;
						//break;
					}
					else
					{
						weight += unit.weight[cand_idx];
					}
				}
                ++j;
                ++k;
            }
            else if(cand_idx < cand_idx2)
                ++j;
            else
                ++k;
        }
        //cout<<"[After] Case 3. return val: "<<weight<<endl;
        return weight;
    }
}

inline long long compute_candidate_size_weight(int u, int num_mapped_neighbors, int un)
{
    bool zero_weight = true;
    int neighbor, neighbor_cand_index, neighbor2, neighbor_cand_index2;
    NodeIndexUnit & unit = indexSet[u];
    if(num_mapped_neighbors == 1){
        temp_intersection_array_index[u][0] = 0;
		neighbor = un;
		neighbor_cand_index = self_pos[neighbor];
		int un_index = query_neighbor_index[u][neighbor];
		for(int j = 0; j < unit.backtrack_size_of_index[un_index][neighbor_cand_index]; ++j){
			//cout<<un_index<<" "<<neighbor_cand_index<<" "<<j<<endl;
			int cand_idx = unit.backtrack_index[un_index][neighbor_cand_index][j];
			if(use_failing_set)
			{
				temp_intersection_array[u][0][temp_intersection_array_index[u][0]++] = cand_idx;
				if(zero_weight && !mapping_flag_data[unit.candidates[cand_idx]]  ) zero_weight = false;
			}
			else
			{
				// if(!mapping_flag_data[unit.candidates[cand_idx]]){
				// 	temp_intersection_array[u][0][temp_intersection_array_index[u][0]++] = cand_idx;
				// }
				temp_intersection_array[u][0][temp_intersection_array_index[u][0]++] = cand_idx;
			}
		}
        if(use_failing_set && zero_weight) return 0;
		else return unit.backtrack_size_of_index[un_index][neighbor_cand_index];
        // else return temp_intersection_array_index[u][0];
    }
    else{
        temp_intersection_array_index[u][num_mapped_neighbors-1] = 0;
        neighbor = un;
        neighbor_cand_index = self_pos[neighbor];
        int un_index = query_neighbor_index[u][neighbor];
        int j = 0, k = 0, cand_idx, cand_idx2, backtrack_size, backtrack_size2;
        backtrack_size = temp_intersection_array_index[u][num_mapped_neighbors-2]; // candidate which is computed previsouly
        backtrack_size2 = unit.backtrack_size_of_index[un_index][neighbor_cand_index]; // candidate which is computed from new query node

        // compute the intersection of two candidates
        while(j < backtrack_size and k < backtrack_size2){
            cand_idx = temp_intersection_array[u][num_mapped_neighbors-2][j];
            cand_idx2 = unit.backtrack_index[un_index][neighbor_cand_index][k];
            if(cand_idx == cand_idx2){
                if(use_failing_set)
                {
        	        temp_intersection_array[u][num_mapped_neighbors-1][temp_intersection_array_index[u][num_mapped_neighbors-1]++] = cand_idx;
                    if(zero_weight and !mapping_flag_data[unit.candidates[cand_idx]]  ) zero_weight = false;
                }
                else
                {
                    if(!mapping_flag_data[unit.candidates[cand_idx]]){
        	            temp_intersection_array[u][num_mapped_neighbors-1][temp_intersection_array_index[u][num_mapped_neighbors-1]++] = cand_idx;
                    }
                }
                ++j;
                ++k;
            }
            else if(cand_idx < cand_idx2)
                ++j;
            else
                ++k;
        }
        if(use_failing_set and zero_weight) return 0;
        else return temp_intersection_array_index[u][num_mapped_neighbors-1];
    }
}


inline double backtrack() {

	long long weight;
	bool flag, all_parents_mapped = true, skip_function_call = false;
	int num_recent_insert[MAX_QUERY_NODE];
	int curr_mark_matching[MAX_QUERY_NODE];
	double count_all_mapping = 0;
	int current, parent_index, pos, data_id, child_data_id, neighbor_data_id, query_id, vertex_in_u;
	int current_overlap=0;
	PointerDagSearchUnit * temp_su;
	NodeIndexUnit * index_unit, * child_index_unit, * neighbor_index_unit;
	fill(num_recent_insert, num_recent_insert + MAX_QUERY_NODE, 0);
	memset(weight_array, 1, sizeof(weight_type) * cnt_node_query);
	memset(mapped_query, false, sizeof(bool) * cnt_node_query);
	memset(mapped_parent, 0, sizeof(int) * cnt_node_query); 
	memset(mapped_neighbors, 0, sizeof(int) * cnt_node_query);
	memset(frontier_index, -1, sizeof(int) * cnt_node_query);
	clear_frontier_improved();
	if(use_failing_set)
	{
		for(int i = 0; i < cnt_node_query; i++)
		{
			adaptive_ancestor[0][i] = bit_vector(cnt_node_query, 0);
		}
	}

	// select backtrack start vertex candidate size of which is smallest.
	// int backtrack_root_node_id = root_node_id;
	int backtrack_root_node_id = -1;
	for(int i = 0; i < cnt_node_query; i++)
	{
		if(NEC_map[i] != -1) continue;
		if(backtrack_root_node_id == -1) backtrack_root_node_id = i;
		else if(indexSet[backtrack_root_node_id].size > indexSet[i].size)
		{
			backtrack_root_node_id = i;
		}	
	}
	// int backtrack_root_node_id = global_matching_order[0];

	if(use_failing_set)
	{
		for(int i = 0; i < cnt_node_query; ++i){
			memset(exist_u[i], false, sizeof(bool) * cnt_node_query);
		}
		memset(ancestor_set_index, 0, sizeof(int) * cnt_node_query);
	}
	mapped_query[backtrack_root_node_id] = true;
	for(int i = query_nodes_array_info[backtrack_root_node_id]; i < query_nodes_array_info[backtrack_root_node_id + 1]; i++)
	{
		int un = query_nodes_array[i];
		if(mapped_query[un]) continue;
		++mapped_neighbors[un];
		if(use_failing_set)
		{
			adaptive_ancestor[mapped_neighbors[un]][un] = adaptive_ancestor[mapped_neighbors[un]-1][un];
			adaptive_ancestor[mapped_neighbors[un]][un] |= adaptive_ancestor[mapped_neighbors[backtrack_root_node_id]][backtrack_root_node_id];
			adaptive_ancestor[mapped_neighbors[un]][un][backtrack_root_node_id] = 1;
		}
	}
	//    cout<<"Candidate size of Root: "<<indexSet[backtrack_root_node_id].size<<endl;

#ifdef MAPPING_FUNCTION_LOG
	cout << "[0] candidates of the root: ";
	for(int i=0;  i<indexSet[backtrack_root_node_id].size; ++i){
		cout<< indexSet[backtrack_root_node_id].candidates[i]<<"("<<i<<"), ";
	}
	cout<<endl;
#endif
	bool found_embedding = false;
	curr_mark_matching[0] = 0;
	memset(indexSet[backtrack_root_node_id].vis, false, sizeof(bool) * max_label_counter);

	for (int i = 0; i < 2 * indexSet[backtrack_root_node_id].size; i++) { // B.1. for each data vertex v in r.C
#ifdef MAPPING_FUNCTION_LOG
		if( current_overlap != 0)
			cerr<<"overlap is not ZERO. Current overlap: "<<current_overlap<<endl;
#endif
		found_embedding = false;
		current_overlap = 0;
		if(b_over_time_limit)
		{
			return count_all_mapping;
		}

		int root_cand_id;
		int idx ;
		if ( i < indexSet[backtrack_root_node_id].size ){
			// first scan
#ifdef FIRSTASC
			idx = i;
#else
			idx = indexSet[backtrack_root_node_id].size - 1 - i;
#endif
			root_cand_id = indexSet[backtrack_root_node_id].candidates[idx];
			if (root_cand_id == -1 || chosen_vertex[root_cand_id])
				continue;
			indexSet[backtrack_root_node_id].vis[idx] = true;
		}
		else{
			// second scan
#ifdef SECONDASC
			idx = i - indexSet[backtrack_root_node_id].size;
#else 
			idx = indexSet[backtrack_root_node_id].size * 2 - 1 - i;
#endif
			root_cand_id = indexSet[backtrack_root_node_id].candidates[idx];
			if (root_cand_id == -1 || indexSet[backtrack_root_node_id].vis[idx])
				continue;

			if(overlap_bound == 0)
				continue;
			else{
				current_overlap++;
			}
		}

		/*
		int root_cand_id = indexSet[backtrack_root_node_id].candidates[i];
		if (root_cand_id == -1)
			continue;

		if(chosen_vertex[root_cand_id]){
			if(overlap_bound == 0)
				continue;
			else{
				current_overlap++;
			}
		}
		*/

		++recursive_call_count;
#ifdef MAPPING_FUNCTION_LOG
		cout<<"rcc: "<<recursive_call_count<<endl; 
		cout << "[0] B.1. for a data vertex v" << root_cand_id << " in r.C(size:"<<indexSet[backtrack_root_node_id].size<<")"<< endl;
#endif


		mapping_flag_data[root_cand_id] = 1;    // B.2. Mark v as visited
		if(use_failing_set)
		{
			mapping_flag_data_int[ root_cand_id ] = backtrack_root_node_id;
		}
		actual_mapping[0] = root_cand_id;       // B.2. M(r) = v
		self_pos[backtrack_root_node_id] = idx;
		char back_trace = 0;
#ifdef MAPPING_FUNCTION_LOG
		cout<<"[0] B.2. M[u"<<backtrack_root_node_id<<"]: v"<<root_cand_id<<", Mark v"<<root_cand_id<<" as visited"<<endl;
#endif
		//must reset this to 1 here
		int current_query_index = 1; // the current query index of the query sequence, because 0 is the root has already been matched
		//mapped[backtrack_root_node_id] = true;

		// B.3. for each neighbor un of r in query graph
		for(int j = query_nodes_array_info[backtrack_root_node_id]; j < query_nodes_array_info[backtrack_root_node_id + 1]; j++)
		{

			int rn = query_nodes_array[j];
			if(NEC_map[rn] != -1) continue;
			if(mapped_query[rn]) continue; // already mapped

#ifdef MAPPING_FUNCTION_LOG
			cout << "["<<current_query_index-1<<"] B.3. for a child of r u" << rn << " in dag(q)"<< endl;
#endif

			if(use_candidate_size_order) weight = compute_candidate_size_weight(rn, mapped_neighbors[rn], backtrack_root_node_id);
			if(use_path_size_order) weight = compute_path_size_weight(rn, mapped_neighbors[rn], backtrack_root_node_id);
			// weight = indexSet[rn].size;


			// become extendable vertex
			if(mapped_neighbors[rn] == 1){
				// if(mapped_neighbors[rn] == DAG_parent_query_size[rn]){
#ifdef MAPPING_FUNCTION_LOG
				cout<<"["<<current_query_index-1<<"] B.4. if all u"<<rn<<"'s parents are mapped then compute_weight("<<rn<<") = "<<weight<<endl;
#endif
				insert_to_frontier_first_improved(rn, weight);
				// insert_to_frontier_first_improved(rn, indexSet[rn].size);
#ifdef MAPPING_FUNCTION_LOG
				cout<<"["<<current_query_index-1<<"] B.5. frontier.insert(u"<<rn<<", weight = "<<weight<<")"<<endl;
#endif
			}
			} // end B.3.

			while (true) 
			{
#ifdef MAPPING_FUNCTION_LOG
				cout << "["<<current_query_index<<"] num recent insert"<< endl;
				for(int i=0; i<current_query_index; ++i){
					cout<<num_recent_insert[i]<<" ";
				}
				cout<<endl;
#endif

				if(cur_k >= num_k || b_over_time_limit)
				{
					ptr_dag_su[current_query_index].address = NULL;
					while (current_query_index != 0)
					{
						current_query_index --;
						mapping_flag_data[ actual_mapping[current_query_index] ] = 0;
						if(use_failing_set)
						{
							mapping_flag_data_int[ actual_mapping[current_query_index] ] = -1;

							while(ancestor_set_index[current_query_index] != 0)
							{
								vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
								exist_u[current_query_index][vertex_in_u] = false; 
							}
						}
						ptr_dag_su[current_query_index].address = NULL; 
					}
					mapping_flag_data[root_cand_id] = 0;
					if(use_failing_set)
						mapping_flag_data_int[root_cand_id] = -1;
#ifdef MAPPING_FUNCTION_LOG
					cout << "["<<current_query_index<<"] A. output M as an embedding of q in G"<< endl;
#endif
					return count_all_mapping; // A. output M as an embedding of q in G
				}// if (b_over_time_limit)



				if (current_query_index == 0)//"No MATCH found!"
					break;
				if (current_query_index == bfs_sequence_index) { // 20170414. found a mapping
					/*
					   assert(cnt_node_query == bfs_sequence_index);
					   if(cnt_node_query != bfs_sequence_index){
					   cerr<<"NOT FULL"<<endl;
					   exit(1);
					   }
					   */

					//======= HERE, we continue to find the full mapping for this core mapping
					// do not return now continue to the roll back process
#ifdef MAPPING_FUNCTION_LOG
					cout<<"found a mapping: exploreCRSequence_indx = "<<exploreCRSequence_indx<<endl;
					cout<<"residual_tree_match_seq_index: "<<residual_tree_match_seq_index<<endl;
					cout<<"NEC_mapping_pair_index: "<<NEC_mapping_pair_index<<endl;
					int temp_cnt = count_all_mapping;
#endif
					if (NEC_mapping_pair_index != 0) // if there are leaves in a query graph
					{
						cout<<"ERROR"<<endl;
						exit(1); //XXX
						//	count_all_mapping += matchLeaves(count_all_mapping); // find mappings of leaf query vertices
						int mlcnt = matchLeaves(count_all_mapping);
						cout<<mlcnt<<endl;
						count_all_mapping += mlcnt; // find mappings of leaf query vertices
					}
					else{
						found_embedding = true;
						count_all_mapping++;
						emb_element& elt = stored_emb[cur_k++];
						for( int idx = 0; idx < bfs_sequence_index; ++idx){
							/*
							   if (chosen_vertex[actual_mapping[idx]] )
							   cerr<<"OVERLAPPED"<<endl;
							   */
							chosen_vertex[actual_mapping[idx]] = true;
							elt.mapped_vertex[idx] = actual_mapping[idx];
						}
					}
#ifdef PRINT_EMB
					if (NEC_mapping_pair_index != 0){ // if there are leaves in a query graph
						int mlcnt=0;
						if(mlcnt = matchLeaves(count_all_mapping)!=0){
							cout<<"matching Leaves: "<<mlcnt<<endl;
							for(int depth_i=0; depth_i < current_query_index; ++depth_i){
								cout<<(int)actual_mapping[depth_i]<<" ";
							}
							cout<<endl;
							//                        cout<<count_all_mapping<<endl;
							assert(0);
						}
					}
					else{
						for(int depth_i=0; depth_i < current_query_index; ++depth_i){
							cout<<(int)actual_mapping[depth_i]<<" ";
						}
						cout<<endl;
						//                        cout<<count_all_mapping<<endl;
					}
#endif
/*
						for(int depth_i=0; depth_i < current_query_index; ++depth_i){
							cout<<(int)actual_mapping[depth_i]<<" ";
						}
						cout<<endl;
*/

#ifdef MAPPING_FUNCTION_LOG
					cout<<"count_all_mapping: "<<count_all_mapping << ", currently found mapping: "<< count_all_mapping - temp_cnt <<endl;
#endif
					if (cur_k >= num_k || (LIMIT > 0 && count_all_mapping >= LIMIT))
					{
						//need to clean up the two array: mapping_flag_data, ptr_dag_su[i].address for next time using
						while (current_query_index != 0)
						{
							current_query_index --;
							mapping_flag_data[ actual_mapping[current_query_index] ] = 0;
							if(use_failing_set)
							{
								mapping_flag_data_int[ actual_mapping[current_query_index] ] = -1;
								while(ancestor_set_index[current_query_index] != 0)
								{
									vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
									exist_u[current_query_index][vertex_in_u] = false; 
								}
							}
							ptr_dag_su[current_query_index].address = NULL; // 20170414
						}
						mapping_flag_data[root_cand_id] = 0;
						if(use_failing_set)
							mapping_flag_data_int[root_cand_id] = -1;
#ifdef MAPPING_FUNCTION_LOG
						cout << "["<<current_query_index<<"] A. output M as an embedding of q in G"<< endl;
#endif
						return count_all_mapping; // A. output M as an embedding of q in G
					} // if (LIMIT > 0 && count_all_mapping >= LIMIT)
					if(use_failing_set)
						mapping_flag_data_int[ actual_mapping[current_query_index - 1] ] = -1;	
					//C.7. Remove just inserted frontier elements
					remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index - 1]); 
					num_recent_insert[current_query_index-1] = 0;
					mapping_flag_data[ actual_mapping[current_query_index - 1] ] = 0; // C.7. Mark v as unvisited
#ifdef MAPPING_FUNCTION_LOG
					cout<<"["<<current_query_index-1<<"] C.7. Mark v"<<actual_mapping[current_query_index - 1]<<" as unvisited."<<endl;
#endif
					current_query_index --;
					if(use_failing_set)
					{
#ifdef MAPPING_FUNCTION_LOG
						int temp_u = ptr_dag_su[current_query_index].vertex;
						cout << "A_prime: "<<one_vector<<", A_prime["<<temp_u<<"]: 1"<<endl;
#endif
						ptr_dag_su[current_query_index].return_vector = one_vector; //20170517
					}
				if(found_embedding){
					if( overlap_bound == cnt_node_query)  continue;

#ifdef PRINT_EMB
					/*
					   for( int idx = 0; idx < bfs_sequence_index; ++idx){
					   cout.width(6);
					   cout<<actual_mapping[idx]<<" ";
					   }
					   cout<<endl;
					   */
#endif
					//TODO: backtrack 위치 조정 for case 1
					int recent_overlap_index = overlap_bound -1;
					
					//while (current_query_index != (recent_overlap_index + 1))
					while ( current_query_index != (recent_overlap_index + 1))
					{
						temp_su = &ptr_dag_su[current_query_index];
						current = temp_su->vertex;
						//
						mapping_flag_data[ actual_mapping[current_query_index] ] = 0;
						if(use_failing_set)
						{
							mapping_flag_data_int[ actual_mapping[current_query_index] ] = -1;

							while(ancestor_set_index[current_query_index] != 0)
							{
								vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
								exist_u[current_query_index][vertex_in_u] = false; 
							}
						}
						ptr_dag_su[current_query_index].address = NULL; 

						//mapped_query

						mapped_query[current] = false;
						for(int i = query_nodes_array_info[current]; i < query_nodes_array_info[current + 1]; i++)
						{
							int un = query_nodes_array[i];
							if(mapped_query[un]) continue;
							--mapped_neighbors[un];
						}

						//queue update
						//remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index]); 
						if(  overlap_bound != 0 ){
#ifdef MAPPING_FUNCTION_LOG
							cout<<"cur "<<cur_k<<endl;
							cout<<"overlap "<<overlap_bound<<endl;
							cout<<"curr "<<current_query_index<<endl;
							cout<<"HEAP: ";
							for(int i=0; i<frontier_node_idx; ++i){
								cout<<frontier_node[i]<<" ";
							}
							cout<<endl;
							cout<<"current "<<current;
							cout<<"num recent insert "<<num_recent_insert[current_query_index]<<endl;
#endif
							reinsert_to_frontier_improved(current, cur_weight[current_query_index], position[current_query_index]);
							remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index-1]);


						}
						num_recent_insert[current_query_index-1] = 0;


						mapping_flag_data[ actual_mapping[current_query_index] ] = 0;
						current_query_index --;

					}
					mapping_flag_data[ actual_mapping[current_query_index] ] = 0;

					found_embedding = false;

					if(   overlap_bound == 0 ){
						//current_overlap--; 
						break;
					}
					else{
						continue;

					}
					//current_overlap = 0;

				}
					continue;
				} // if (current_query_index == bfs_sequence_index)

				temp_su = &ptr_dag_su[current_query_index];
				// Find \intersection_{u_p in u.p} (N^{u_p}_{u}(M(u_p)))
				if (temp_su->address == NULL) { 
					++recursive_call_count;
#ifdef MAPPING_FUNCTION_LOG
					cout<<"rcc!: "<<recursive_call_count<<endl; 
#endif
					// current = global_matching_order[current_query_index];
					pop_from_frontier_improved(current, cur_weight[current_query_index], position[current_query_index]);
#ifdef MAPPING_FUNCTION_LOG
					cout << "["<<current_query_index<<"] C.1. Pop a vertex with the minimum weight: " << current << endl;
#endif
					temp_su->vertex = current;
					index_unit = &indexSet[current];
					mapped_query[current] = true; 
					for(int i = query_nodes_array_info[current]; i < query_nodes_array_info[current + 1]; i++)
					{
						int un = query_nodes_array[i];
						if(mapped_query[un]) continue;
						++mapped_neighbors[un];
						if(use_failing_set)
						{
							adaptive_ancestor[mapped_neighbors[un]][un] = adaptive_ancestor[mapped_neighbors[un]-1][un];
							adaptive_ancestor[mapped_neighbors[un]][un] |= adaptive_ancestor[mapped_neighbors[current]][current];
							adaptive_ancestor[mapped_neighbors[un]][un][current] = 1;
						}
					}

					if(use_failing_set)
						temp_su->return_vector = zero_vector;

					temp_su->address = temp_intersection_array[current][mapped_neighbors[current]-1];
					temp_su->address_size = temp_intersection_array_index[current][mapped_neighbors[current]-1];
#ifdef MAPPING_FUNCTION_LOG
					cout<<mapped_neighbors[current]<<endl;
					string temp ="";
					for(int i = 0; i < temp_su->address_size; ++i){
						temp += to_string(index_unit->candidates[temp_su->address[i]])+"("+to_string(temp_su->address[i])+"), ";
					}
					cout << "["<<current_query_index<<"] C.2. Find intersection of of N^{u.p}_{u}(M(u_p)) for all u_p in u.p: " << temp <<"(size: "<<temp_su->address_size<<")"<< endl;
					cout << "read skip_function_call: " << (skip_function_call?true:false) << endl;
#endif
					if (temp_su->address_size == 0 || skip_function_call)
					{
						temp_su->address = NULL;
						if(use_failing_set)
							temp_su->firstly_visited = true;
						current_query_index--; // roll back one node in the matching sequence
						skip_function_call = false;


						if (current_query_index != 0)
						{
							// C.8. frontier.insert(u, weight)
							reinsert_to_frontier_improved(current, cur_weight[current_query_index+1], position[current_query_index+1]); 
#ifdef MAPPING_FUNCTION_LOG
							cout <<"["<<current_query_index+1<<"] C.8. frontier.insert("<<current<<")!"<<endl;
#endif

							if(use_failing_set)
							{
								if(!skip_function_call){
									while(ancestor_set_index[current_query_index+1] != 0){
										vertex_in_u = ancestor_set[current_query_index+1][--ancestor_set_index[current_query_index+1]];
										temp_su->return_vector |= adaptive_ancestor[mapped_neighbors[vertex_in_u]][vertex_in_u];
										exist_u[current_query_index+1][vertex_in_u] = false; 
									}
								}
#ifdef MAPPING_FUNCTION_LOG
								cout << "\tReturn: "<<temp_su->return_vector<<endl;
#endif
							}
						} // if (current_query_index != 0)

						mapped_query[current] = false;
						for(int i = query_nodes_array_info[current]; i < query_nodes_array_info[current + 1]; i++)
						{
							int un = query_nodes_array[i];
							if(mapped_query[un]) continue;
							--mapped_neighbors[un];
						}
						continue;
					} // if (temp_su->address_size == 0 || skip_function_call)
					temp_su->address_pos = 0;
					while(temp_su->address_pos < 2* temp_su->address_size){
						if (temp_su->address_pos < temp_su->address_size){
#ifdef FIRSTASC
							int idx = temp_su->address_pos;
#else
							int idx = temp_su->address_size - 1 - temp_su->address_pos;
#endif
							int pos = temp_su->address[idx];
							int v = index_unit->candidates[pos];
							if( !chosen_vertex[v])
								break;
						}
						else{
#ifdef SECONDASC
							int idx = temp_su->address_pos - temp_su->address_size;
#else
							int idx = 2 * temp_su->address_size - 1 - temp_su->address_pos;
#endif
							int pos = temp_su->address[idx];
							int v = index_unit->candidates[pos];
							if( chosen_vertex[v])
								break;
						}
						++temp_su->address_pos;
					}
				} 
				else // temp_su->address != NULL
				{ 
					current = temp_su->vertex;
					index_unit = &indexSet[current];
					temp_su->address_pos++; // update the index by one

					while(temp_su->address_pos < 2* temp_su->address_size){
						if (temp_su->address_pos < temp_su->address_size){
#ifdef FIRSTASC
							int idx = temp_su->address_pos;
#else
							int idx = temp_su->address_size - 1 - temp_su->address_pos;
#endif
							int pos = temp_su->address[idx];
							int v = index_unit->candidates[pos];
							if( !chosen_vertex[v])
								break;
						}
						else{
#ifdef SECONDASC
							int idx = temp_su->address_pos - temp_su->address_size;
#else
							int idx = temp_su->address_size * 2 - 1 - temp_su->address_pos ;
#endif
							int pos = temp_su->address[idx];
							int v = index_unit->candidates[pos];
							if( chosen_vertex[v])
								break;
						}
						++temp_su->address_pos;
					}

					if ( temp_su->address_pos == 2 * temp_su->address_size || skip_function_call) 
					{
						// C.8. frontier.insert(u, weight). weight_array[current] does not need to be changed
						reinsert_to_frontier_improved(current,cur_weight[current_query_index],position[current_query_index]); 
#ifdef MAPPING_FUNCTION_LOG
						cout<<"["<<current_query_index<<"] C.8. frontier.insert("<<current<<")!!"<<endl; 
						cout << "read skip_function_call: " << (skip_function_call?true:false) << endl;
#endif
						if(use_failing_set)
						{
							if(!skip_function_call)
							{
								while(ancestor_set_index[current_query_index] != 0)
								{
									vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
									temp_su->return_vector |= adaptive_ancestor[mapped_neighbors[vertex_in_u]][vertex_in_u];
									exist_u[current_query_index][vertex_in_u] = false; 
								}
							}
#ifdef MAPPING_FUNCTION_LOG
							cout << "\tReturn: "<<temp_su->return_vector<<endl;
#endif
						}
						mapped_query[current] = false;
						for(int i = query_nodes_array_info[current]; i < query_nodes_array_info[current + 1]; i++)
						{
							int un = query_nodes_array[i];
							if(mapped_query[un]) continue;
							--mapped_neighbors[un];
						}
						if(use_failing_set)
						{
							skip_function_call = false;
							temp_su->firstly_visited = true;
						}
						temp_su->address = NULL;
						current_query_index--;
						if(use_failing_set)
							mapping_flag_data_int[ actual_mapping[ current_query_index ] ] = -1;
						//C.7. Remove just inserted frontier elements
						remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index]); 
						num_recent_insert[current_query_index] = 0;
						mapping_flag_data[ actual_mapping[ current_query_index ] ] = 0;
						if( chosen_vertex[ actual_mapping[ current_query_index ] ])
							--current_overlap;

#ifdef MAPPING_FUNCTION_LOG
						cout<<"["<<current_query_index<<"] C.7. Mark v"<<actual_mapping[current_query_index]<<" as unvisited, current_query_index: "<<endl;
#endif
						if(use_failing_set)
						{
							bit_vector& A_prime = ptr_dag_su[current_query_index+1].return_vector;
#ifdef MAPPING_FUNCTION_LOG
							int temp_u = ptr_dag_su[current_query_index].vertex;
							cout << "A_prime: "<<A_prime<<", A_prime["<<temp_u<<"]: "<<A_prime[temp_u]<<endl;
#endif
							if(current_query_index != 0)
							{
								if(A_prime[ptr_dag_su[current_query_index].vertex] == 0)
								{
									ptr_dag_su[current_query_index].return_vector = A_prime;
									skip_function_call = true; //20170517
								}
								else
								{
									ptr_dag_su[current_query_index].return_vector |= A_prime;
								}
							}
#ifdef MAPPING_FUNCTION_LOG
							cout << "write skip_function_call: " << (skip_function_call?true:false) << endl;
							cout << "A: "<<ptr_dag_su[current_query_index].return_vector<<endl;
#endif
						}
						continue;
					} // if ( temp_su->address_pos == temp_su->address_size || skip_function_call) 
				} // if (temp_su->address == NULL) else

				back_trace = 0; //actually, this line is not necessary, when processed here, the back_trace must be false...
				while (true) {
					//break, until find a mapping for this node
					//or cannot find a mapping after examining all candidates
					//no recursive call count increase here!
					if(b_over_time_limit)
					{
						if(use_failing_set)
						{
							while(ancestor_set_index[current_query_index] != 0)
							{
								vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
								exist_u[current_query_index][vertex_in_u] = false; 
							}
						}
						ptr_dag_su[current_query_index].address = NULL; // 20170414
						while (current_query_index != 0)
						{
							current_query_index --;
							mapping_flag_data[ actual_mapping[current_query_index] ] = 0;
							if(use_failing_set)
							{
								mapping_flag_data_int[ actual_mapping[current_query_index] ] = -1;
								while(ancestor_set_index[current_query_index] != 0)
								{
									vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
									exist_u[current_query_index][vertex_in_u] = false; 
								}
							}	
							ptr_dag_su[current_query_index].address = NULL; // 20170414
						}
						mapping_flag_data[root_cand_id] = 0;
						if(use_failing_set)
							mapping_flag_data_int[root_cand_id] = -1;   
#ifdef MAPPING_FUNCTION_LOG
						cout <<"["<<current_query_index<<"] A. output M as an embedding of q in G"<< endl;
#endif
						return count_all_mapping; // A. output M as an embedding of q in G
					} // if (b_over_time_limit)

					if(temp_su->address_pos < temp_su->address_size){
#ifdef FIRSTASC
						pos = temp_su->address[temp_su->address_pos]; // the position of v in C_M(u)
#else
						pos = temp_su->address[temp_su->address_size - 1 - temp_su->address_pos]; // the position of v in C_M(u)
#endif
					}
					else{
#ifdef SECONDASC
						pos = temp_su->address[temp_su->address_pos - temp_su->address_size]; // the position of v in C_M(u)
#else
						pos = temp_su->address[temp_su->address_size * 2 - 1 - temp_su->address_pos ]; // the position of v in C_M(u)
#endif
					}
					data_id = index_unit->candidates[pos]; // v\in C_M(u)
#ifdef MAPPING_FUNCTION_LOG
					cout<<"TEMP: temp_su->vertex: "<<temp_su->vertex<<", temp_su->address_pos: "<<temp_su->address_pos<<", temp_su->address_size: "<<temp_su->address_size<<", data_id: v"<<data_id<<endl;
#endif
					if (!mapping_flag_data[data_id]) // C.3. if v is unvisited 
					{ 
						if(chosen_vertex[data_id]){
#ifdef MAPPING_FUNCTION_LOG
							cout<<"["<<current_query_index<<"] OV.1. if data vertex v"<<data_id<<" is already chosen"<<endl;
							cout<<"number recent "<<num_recent_insert[current_query_index]<<endl;
#endif
							if(  current_overlap >= overlap_bound){
#ifdef MAPPING_FUNCTION_LOG
								cout<<"["<<current_query_index<<"] OV.2. fully overlapped"<<endl;
#endif
								temp_su->address_pos++;
								if (temp_su->address_pos == 2 * temp_su->address_size) 
								{ // no more data id, so cannot find a match for this query node
									back_trace = 1; // indicate that no result is being found, so we need to trace back_trace
									break;
								}
								else
									continue;
							}
							else{
								++current_overlap;
							}
						}

#ifdef MAPPING_FUNCTION_LOG
						cout<<"["<<current_query_index<<"] C.3. if data vertex v"<<data_id<<" is unvisited"<<endl;
#endif
						actual_mapping[current_query_index] = data_id; // C.4. M(u) = v
						self_pos[current] = pos;
						mapping_flag_data[data_id] = 1; // C.4. Mark v as visited

						if(use_failing_set)
							mapping_flag_data_int[data_id] = current; // C.4. Mark v as visited

#ifdef MAPPING_FUNCTION_LOG
						cout<<"["<<current_query_index<<"] C.4. M[u"<<current<<"]: v"<<data_id<<", Mark v"<<data_id<<" as visited"<<endl;
#endif
						if(use_path_size_order)
						{
							frontier_min_index = -1;
							frontier_min_weight = WEIGHT_MAX;
						}

						// C.5. for each neighbor un of u in query graph s.t. un is extendable vertex
						flag = true;
						for(int j = query_nodes_array_info[current]; j < query_nodes_array_info[current + 1]; j++)
						{
							int current_neighbor = query_nodes_array[j];
							if(NEC_map[current_neighbor] != -1) continue;
							if(mapped_query[current_neighbor]) continue;

							if(use_candidate_size_order) weight = compute_candidate_size_weight(current_neighbor, mapped_neighbors[current_neighbor], current);
							if(use_path_size_order) weight = compute_path_size_weight(current_neighbor, mapped_neighbors[current_neighbor], current);


							// if(DAG_parent_query_size[current_neighbor] > 1 && mapped_neighbors[current_neighbor] == 1) continue;
#ifdef MAPPING_FUNCTION_LOG
							cout << "["<<current_query_index<<"] C.5. for u"<<current<<"'s child u"<<current_neighbor<<" in dag(q) s.t. all parents were mapped, compute_weight("<<current_neighbor<<") = "<<weight<<endl;
#endif


							if(weight == 0)
							{
								flag = false;

								if(use_failing_set)
								{
									if(!exist_u[current_query_index][current_neighbor])
									{
										ancestor_set[current_query_index][ancestor_set_index[current_query_index]++] = current_neighbor;
										exist_u[current_query_index][current_neighbor] = true;
#ifdef MAPPING_FUNCTION_LOG
										cout<<"U = U u {u"<<current_neighbor<<"}!"<<endl;
#endif
									}

									neighbor_index_unit = &indexSet[current_neighbor];
									for(int i = 0; i < temp_intersection_array_index[current_neighbor][mapped_neighbors[current_neighbor]-1]; ++i)
									{
										int neighbor_cand_id = temp_intersection_array[current_neighbor][mapped_neighbors[current_neighbor]-1][i];
										neighbor_data_id = neighbor_index_unit->candidates[neighbor_cand_id];

										if(mapping_flag_data[neighbor_data_id])
										{
											query_id = mapping_flag_data_int[neighbor_data_id];
											temp_su->return_vector[query_id] = 1;

											if(!exist_u[current_query_index][query_id])
											{
												ancestor_set[current_query_index][ancestor_set_index[current_query_index]++] = query_id;
												exist_u[current_query_index][query_id] = true;
											}

#ifdef MAPPING_FUNCTION_LOG
											cout<<"U = U u {u"<<query_id<<"}!!"<<endl;
#endif
										}
									}
								}

#ifdef MAPPING_FUNCTION_LOG
								cout << "["<<current_query_index<<"] C.5f. flag = false"<<endl;
#endif
								break;
							}

							// become extendable vertex
							if(mapped_neighbors[current_neighbor] == 1)
								// if(mapped_neighbors[current_neighbor] == DAG_parent_query_size[current_neighbor])
							{
								insert_to_frontier_first_improved(current_neighbor, weight);
								// insert_to_frontier_first_improved(current_neighbor, indexSet[current_neighbor].size);
								num_recent_insert[current_query_index]++;

#ifdef MAPPING_FUNCTION_LOG
								cout<<"["<<current_query_index<<"] C.6(1). frontier.insert(u"<<current_neighbor<<", weight = "<<weight<<"), current_query_index: "<<current_query_index<<endl;
#endif
							}
							// update weight if current_neighbor is an extendable vertex
							else if(0 <= frontier_index[current_neighbor] && frontier_index[current_neighbor] < frontier_node_idx)
							{
								weight_array[frontier_index[current_neighbor]] = weight;
								// weight_array[frontier_index[current_neighbor]] = indexSet[current_neighbor].size;
#ifdef MAPPING_FUNCTION_LOG
								cout<<"["<<current_query_index<<"] C.6(2). frontier.update weight(u"<<current_neighbor<<", weight = "<<weight<<"), current_query_index: "<<current_query_index<<endl;
#endif
							}
						} // end C.5.

						if(!flag)
						{
							temp_su->address_pos++;
							// C.7. Remove just inserted frontier elements
							remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index]); 
							num_recent_insert[current_query_index] = 0;
							mapping_flag_data[ actual_mapping[ current_query_index ] ] = 0; // C.7. Mark v as unvisited 
							if( chosen_vertex[ actual_mapping[ current_query_index ] ])
								--current_overlap;

							if(use_failing_set)
								mapping_flag_data_int[ actual_mapping[ current_query_index ] ] = -1;
#ifdef MAPPING_FUNCTION_LOG
							cout<<"["<<current_query_index<<"] C.7f. Mark v"<<actual_mapping[current_query_index]<<" as unvisited."<<endl;
#endif
							flag = true;
							if (temp_su->address_pos == 2 * temp_su->address_size)
							{
								back_trace = 1; 
								break;
							}
							else
								continue;
						}
						break;
					} 
					else //mapping NOT OK! // if (mapping_flag_data[data_id]) 
					{ 
						if(use_failing_set)
						{
							if(!exist_u[current_query_index][current])
							{
								ancestor_set[current_query_index][ancestor_set_index[current_query_index]++] = current;
								exist_u[current_query_index][current] = true;      
#ifdef MAPPING_FUNCTION_LOG
								cout<<"U = U u {u"<<current<<"}!!!!"<<endl;
#endif
							} 
							query_id = mapping_flag_data_int[data_id];
							temp_su->return_vector[query_id] = 1;
							if(!exist_u[current_query_index][query_id])
							{
								ancestor_set[current_query_index][ancestor_set_index[current_query_index]++] = query_id;
								exist_u[current_query_index][query_id] = true;      
#ifdef MAPPING_FUNCTION_LOG
								cout<<"U = U u {u"<<query_id<<"}"<<endl;
#endif
							} 
						}
						temp_su->address_pos++;//not ok, then we need the next result
						if (temp_su->address_pos == 2 * temp_su->address_size) 
						{ // no more data id, so cannot find a match for this query node
							back_trace = 1; // indicate that no result is being found, so we need to trace back_trace
							break;
						}
					}
				} //end while
#ifdef MAPPING_FUNCTION_LOG
				cout<<"back_trace: "<<(back_trace?"true":"false")<<", flag: "<<(flag?"true":"false")<<endl;
#endif
				if (back_trace) //Backtrack!
				{ 
					back_trace = 0;
					// C.8. frontier.insert(u, weight). weight_array[current] does not need to be changed
					reinsert_to_frontier_improved(current,cur_weight[current_query_index],position[current_query_index]); 
#ifdef MAPPING_FUNCTION_LOG
					cout<<"["<<current_query_index<<"] C.8. frontier.insert("<<temp_su->vertex<<")!!!"<<endl; 
					cout<<"number recent insert "<<num_recent_insert[current_query_index]<<endl;
#endif
					if(use_failing_set)
					{
						if(!skip_function_call)
						{
							while(ancestor_set_index[current_query_index] != 0)
							{
								vertex_in_u = ancestor_set[current_query_index][--ancestor_set_index[current_query_index]];
								temp_su->return_vector |= adaptive_ancestor[mapped_neighbors[vertex_in_u]][vertex_in_u];
								exist_u[current_query_index][vertex_in_u] = false; 
							}
						}
#ifdef MAPPING_FUNCTION_LOG
						cout << "\tReturn: "<<temp_su->return_vector<<endl;
#endif
					}               
					mapped_query[current] = false;
					for(int i = query_nodes_array_info[current]; i < query_nodes_array_info[current + 1]; i++)
					{
						int un = query_nodes_array[i];
						if(mapped_query[un]) continue;
						--mapped_neighbors[un];
					}
					if(use_failing_set)
					{
						skip_function_call = false;
						temp_su->firstly_visited = true;
					}
					temp_su->address = NULL;
					current_query_index--;
					if(use_failing_set)
						mapping_flag_data_int[ actual_mapping[ current_query_index ] ] = -1;
					// C.7. Mark v as unvisited
					mapping_flag_data[ actual_mapping[ current_query_index ] ] = 0; 
					if( chosen_vertex[ actual_mapping[ current_query_index ] ])
						--current_overlap;
					// C.7. Remove just inserted frontier elements
					remove_recently_inserted_from_frontier_improved(num_recent_insert[current_query_index]); 
					num_recent_insert[current_query_index] = 0;
#ifdef MAPPING_FUNCTION_LOG
					cout<<"["<<current_query_index<<"] C.7. Mark v"<<actual_mapping[current_query_index]<<" as unvisited."<<endl;
#endif
					if(use_failing_set)
					{
						bit_vector& A_prime = ptr_dag_su[current_query_index+1].return_vector;
#ifdef MAPPING_FUNCTION_LOG
						int temp_u = ptr_dag_su[current_query_index].vertex;
						cout << "A_prime: "<<A_prime<<", A_prime["<<temp_u<<"]: "<<A_prime[temp_u]<<endl;
#endif
						if(A_prime[ptr_dag_su[current_query_index].vertex] == 0)
						{
							ptr_dag_su[current_query_index].return_vector = A_prime;
							skip_function_call = true; // 20170517
						}
						else
						{
							ptr_dag_su[current_query_index].return_vector |= A_prime;
						}                          
#ifdef MAPPING_FUNCTION_LOG
						cout << "write skip_function_call: " << (skip_function_call?true:false) << endl;
						cout << "A: "<<ptr_dag_su[current_query_index].return_vector<<endl;
#endif
					}
					// modified
				} 
				else // if (!back_trace)
				{
					current_query_index++;
				}
			} // while
			//===========================================================================================
			mapping_flag_data[root_cand_id] = 0;
			clear_frontier_improved(); // clear frontier
#ifdef MAPPING_FUNCTION_LOG
			cout<<"["<<current_query_index<<"] B.6. Clear frontier. Mark v"<<root_cand_id<<" as unvisited, current_query_index: "<<current_query_index<<endl;
#endif
			if(use_failing_set)
			{
				mapping_flag_data_int[root_cand_id] = -1;
				skip_function_call = false;
			}
		} // for indexSet[root]    
		mapped_query[backtrack_root_node_id] = false;
		for(int i = query_nodes_array_info[backtrack_root_node_id]; i < query_nodes_array_info[backtrack_root_node_id + 1]; i++)
		{
			int un = query_nodes_array[i];
			if(mapped_query[un]) continue;
			--mapped_neighbors[un];
		}

		return count_all_mapping;
	}

inline double backtrack_overall(){
/*
    // find embeddings with no overlapped vertex
    overlap_bound = 0;
    double num_emb_0 = backtrack();
    // find embeddings with one overlapped vertex
    //double num_emb_1 = 0;
    overlap_bound = 1;
    double num_emb_1 = backtrack();

    overlap_bound = 2;
    double num_emb_2 = backtrack();

    overlap_bound = 3;
    double num_emb_3 = backtrack();
//    cout<<num_emb_0<< " "<<num_emb_1<<endl;
    return num_emb_0 + num_emb_1 + num_emb_2 + num_emb_3;
    */
    double num_emb = 0;
    for(overlap_bound = 0; overlap_bound <= cnt_node_query; ++overlap_bound){
//cout<<"-----------"<<overlap_bound<<endl;
        num_emb += backtrack();

        if(cur_k >= num_k || b_over_time_limit || num_emb >= LIMIT)
        break;
    }
/*
    if(cur_k < num_k && !b_over_time_limit && num_emb < LIMIT){
	    phase1 = false;
	    num_emb += backtrack();
    }
*/
    return num_emb;
}

inline void transfer_TurboISOQuery_to_myquery() {

    vector<vector<int> > nodes_query; //the array store a node's adjacent list of query graph
	ifstream fin(querygraphFileFolder);
	string line;
	vector<string> v;
	double degree_query = 0;
	int id = 0;
	int node_index = 0;

	getline(fin, line);

	while (getline(fin, line)) {

		if (line.at(0) == 't') {

			cout << "t " << id << " " << cnt_node_query << " " << degree_query * 2 << endl;
			for (int i = 0; i < cnt_node_query; i++){
				cout << i << " " << nodes_label_query[i] << " " << nodes_query[i].size() << " ";
				for (int j = 0; j < nodes_query[i].size(); j++){
					cout << nodes_query[i][j] << " ";
				}
				cout << endl;
			}

			nodes_query.clear();
			node_index = 0;
			cnt_node_query = 0;
			degree_query = 0;
			id ++;
		}

		if (line.at(0) == 'v') {
			vector <int> node;
			nodes_query.push_back(node);
			split(line, v, ' ');
			int label = atoi(v[2].c_str());
			nodes_label_query[node_index] = label;
			node_index++;
			cnt_node_query ++;
		}

		if (line.at(0) == 'e') {
			split(line, v, ' ');
			int int_node_left = atoi(v[1].c_str());
			int int_node_right = atoi(v[2].c_str());
			//put the nodes into the adjacent list
			nodes_query[int_node_left].push_back(int_node_right);
			nodes_query[int_node_right].push_back(int_node_left);
			degree_query++;
		}

	}
	cout << "t " << id << " " << cnt_node_query << " " << degree_query * 2 << endl;
	for (int i = 0; i < cnt_node_query; i++){
		cout << i << " " << nodes_label_query[i] << " " << nodes_query[i].size() << " ";
		for (int j = 0; j < nodes_query[i].size(); j++)
			cout << nodes_query[i][j] << " ";
		cout << endl;
	}

	fin.close();
}

void print_DAG()
{
	cerr << "-----------DAG child-----------" << endl;
	for(int i = 0; i < cnt_node_query; i++)
	{
		cerr << i << " th node: ";
		if( DAG_child_query_size[i] == 0 )
			cerr << "leaf";

		for(int j = 0; j < DAG_child_query_size[i]; j++)
		{
			cerr << DAG_child_query[i][j];
			if( j != DAG_child_query_size[i] - 1)
				cerr << ", ";
		}
		cerr << endl;
	}

	cerr << "-----------DAG parent-----------" << endl;
	for(int i = 0; i < cnt_node_query; i++)
	{
		cerr << i << " th node: ";
		if( DAG_parent_query_size[i] == 0 )
			cerr << "root";
		for(int j = 0; j < DAG_parent_query_size[i]; j++)
		{
			cerr << DAG_parent_query[i][j];
			if( j != DAG_parent_query_size[i] - 1)
				cerr << ", ";
		}
		cerr << endl;
	}
	cerr << "-----------DAG ancestor---------" << endl;
	for(int i = 0; i < cnt_node_query; i++)
	{
		cerr << i << " th node: ";
		cerr << DAG_ancestor[i] << endl;
	}
	cerr << "-----------DAG print END-----------" << endl << flush;
}

//this function prints lavel and degree for each node in DAG of query graph
void print_DAG_label_and_degree()
{
	//label uses [1-cnt_unique_label]
	/*
	for(int i = 1; i < cnt_unique_label+1; i++)
	{
		int frequency = label_freqency[i];
		int frequency_rank = label_frequency_rank[i];
		cout << i << ": freq=" << frequency << " rank=" << frequency_rank << endl << flush;
	}
	*/
	for(int i = 0; i < cnt_node_query; i++)
	{
		int node = simulation_sequence[i];
		int level = BFS_level_query[node];
		int label = nodes_label_query[node];
		int degree = node_degree_query[node];

		cout << node << ": level=" << level << " label=" << label << " degree=" << degree << endl << flush; 
	}
}

long long total_candidate_count = 0;
long long getCandidateCount()
{
	long long sum = 0;
	for(int i = 0; i < cnt_node_query; i++)
	{
		int current_node = i;

		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node)
		{
			NodeIndexUnit& cur_node_unit = indexSet[NEC_map[current_node]];
			long long size = 0;
			for(int j = 0; j < cur_node_unit.size; j++)
			{
				if(cur_node_unit.candidates[j] != -1)
					size++;
			}
			sum += size;
			continue;
		}

		NodeIndexUnit& cur_node_unit = indexSet[current_node];
		long long size = 0;
		for(int j = 0; j < cur_node_unit.size; j++)
		{
			if(cur_node_unit.candidates[j] != -1)
				size++;
		}
		sum += size;
	}
	return sum;
}

void print_CPI()
{
	cerr << "------------CORE nodes-----------" << endl;
	for(int i = 0; i < cnt_node_query ; i++)
	{
		if( core_number_query[i] > 1 )
			cerr << i << ", ";
	}
	cerr << endl;


	cerr << "------------CPI start------------" << endl;
	cerr << "---------------candidates--------------" << endl;

	for(int i = 0; i < cnt_node_query; i++)
	{
		int current_node = i;

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node)
		{
			cerr << current_node << " th node's candidates: same as u" << NEC_map[current_node] << " (NEC)" << endl;
			continue;
		}

		NodeIndexUnit& cur_node_unit = indexSet[current_node];
		//cerr << current_node << " th node's candidates(size: "<<cur_node_unit.size<<"): ";
		cerr << current_node << " th node's candidates(size: ";
		int size = 0;
		for(int j = 0; j < cur_node_unit.size; j++)
		{
			if(cur_node_unit.candidates[j] != -1)
				size++;
		}
		cerr << size << "): ";
		for(int j = 0; j < cur_node_unit.size; j++)
		{
			if(cur_node_unit.candidates[j] == -1)
				continue;

			cerr << cur_node_unit.candidates[j] << "("<<j<<"), ";
		}
		cerr << endl;
	}


	cerr << "---------------N array--------------" << endl;

    
	for(int i = 0; i < cnt_node_query; i++)
	{
		int current_node = i;

		//NEC boost
		if(NEC_map[current_node] != -1 && NEC_map[current_node] != current_node)
		{
                        int parent_node = DAG_parent_query[current_node][0];
			cerr << "N^{u"<<parent_node<<"}_{u" << current_node << "]: same as N^{u"<<parent_node<<"}_{u" << NEC_map[current_node] << "} (NEC)" <<endl; 
			continue;
		}

		for(int dag_parent_index = 0; dag_parent_index < DAG_parent_query_size[current_node]; dag_parent_index++)
		{
			int parent_node = DAG_parent_query[current_node][dag_parent_index];	
			NodeIndexUnit& cur_node_unit = indexSet[current_node];
			NodeIndexUnit& parent_node_unit = indexSet[parent_node];

			for(int j = 0; j < parent_node_unit.size; j++)
			{
				if(parent_node_unit.candidates[j] == -1)
					continue;

				//cerr << "|N^{u" << parent_node << "}_{u" << current_node << "}|: "<<cur_node_unit.backtrack_size_of_index[0][j]<<endl;; 
				//cerr << "N^{u" << parent_node << "}_{u" << current_node << "} (v" << parent_node_unit.candidates[j] << "): "; 
				//for(int k = 0; k < cur_node_unit.backtrack_size_of_index[dag_parent_index][j]; k++)
				//{
				//	cerr << "v" << cur_node_unit.candidates[ cur_node_unit.backtrack_index[dag_parent_index][j][k] ];
				//	if(k != cur_node_unit.backtrack_size_of_index[dag_parent_index][j] - 1)
				//		cerr << ", ";
				//}
				//cerr << endl;
				//<sort test>
				
				int prev = -100;
				for(int k = 0; k < cur_node_unit.backtrack_size_of_index[dag_parent_index][j]; k++)
				{
					if( prev > cur_node_unit.backtrack_index[dag_parent_index][j][k] )
					{
						cout << "[ERROR] not sorted" << endl; 
					}
					prev = cur_node_unit.backtrack_index[dag_parent_index][j][k];
				}
				cerr << "N^{u" << parent_node << "}_{u" << current_node << "} (v" << parent_node_unit.candidates[j] << "): "; 
				for(int k = 0; k < cur_node_unit.backtrack_size_of_index[dag_parent_index][j]; k++)
				{
					cerr << cur_node_unit.backtrack_index[dag_parent_index][j][k];
					if(k != cur_node_unit.backtrack_size_of_index[dag_parent_index][j] - 1)
						cerr << ", ";
				}
				cerr << endl;
				
                                //</sort test>
			}
		}
	}
    
    
	/*
	cerr << "---------------N array--------------" << endl;

	for(int i = 0; i < cnt_node_query; i++)
	{
		int current_node = i;
		int parent_node = BFS_parent_query[current_node];
		
		NodeIndexUnit& cur_node_unit = indexSet[current_node];
		NodeIndexUnit& parent_node_unit = indexSet[parent_node];

		for(int j = 0; j < parent_node_unit.size; j++)
		{
			cerr << "N[" << current_node << "][" << parent_node << "] (" << j << "): "; 
			for(int k = 0; k < cur_node_unit.size_of_index[ j ]; k++)
			{
				cerr << cur_node_unit.index[j][k] << ", ";
			}
			cerr << endl;
		}
	}*/
	/*
	cerr << "----------------path----------------" << endl;

	for(int i = 0; i < cnt_node_query; i++)
	{
		int current_node = i;
		
		NodeIndexUnit& cur_node_unit = indexSet[current_node];
		cerr << current_node << " th node's path: ";
		for(int j = 0; j < cur_node_unit.size; j++)
		{
			cerr << cur_node_unit.path[j] << ", ";
		}
		cerr << endl;
	}
	*/
	cerr << "------------CPI end--------------" << endl << flush;
}

int computeCoverage(){
		memset(marking_cnt, 0, sizeof(int)*cnt_node);

//cout<<"Embeddings "<<endl;
		for(int i=0; i<cur_k; ++i){
				emb_element& elt = stored_emb[i];
				for(int depth_i=0; depth_i < cnt_node_query; ++depth_i){
//cout<<elt.mapped_vertex[depth_i] <<" ";
						marking_cnt[ elt.mapped_vertex[depth_i] ] ++;
				}
//cout<<endl;
		}

		int res = 0;
		for(int  i =0; i<cnt_node; ++i){
				if(marking_cnt[i] > 0) ++res;
		}
		return res;
}

//Usage 1: run DAF
//./program datafile queryfile #query #matches
//default:  datagraph label <= 100: path-size order
//          datagraph label > 100: candidate-size order
//          use failing set always
//option:   -p: path-size order         //use path-size order
//          -c: candidate-size order    //use candidate-size order
//          -f: disable failing set     //not use failing set
//          -d: datagraph
//          -q: queryfile
//          -m: #match
//          -n: #query
//Usage 2: transfer turboiso_query to DAF_query
//./program -t turboiso_queryfile
int main(int argc, char *argv[]) {

    for(int i = 1; i < argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
                case 't':
		            querygraphFileFolder = argv[i+1];
            		transfer_TurboISOQuery_to_myquery();
                    return 0;
                case 'p':
                    use_path_size_order = true;
                    use_candidate_size_order = false;
                    order_flag = true;
                    break;
                case 'c':
                    use_path_size_order = false;
                    use_candidate_size_order = true;
                    order_flag = true;
                    break;
                case 'f':
                    use_failing_set = false;
                    break;
                case 'd':
	                datagraphFile = argv[i+1];
                    break;
                case 'q':
	                querygraphFileFolder = argv[i+1];
                    break;
                case 'm':
                    getLimit_full(argv[i+1], LIMIT);
                    break;
                case 'n':
                    count_query_file = atoi(argv[i+1]);
                    break;
				case 'k':
					num_k = atoi(argv[i+1]);
					break;
            }
        }
    }
    /*
	cout << "Data File :" << datagraphFile << endl;
	cout << "Query file:" << querygraphFileFolder << endl;
    cout<<"LIMIT: " << LIMIT << endl;
    */

	cerr << "Reading data graph ...";
	preprocessDataGraph(datagraphFile);
	readAndProcessDataGraph(datagraphFile);
	cerr << "OK!" << endl;
	cerr << "Building index...";
	coreDecomposition_data();
	parametersInitilisingBeforeQuery();
	initializeStatisticParameters();
	cerr << "OK!" << endl;

	fin_query.open(querygraphFileFolder);
	int cnt_end_query = 0;

	marking_cnt = new int[cnt_node];

	for (int i = 0; i < count_query_file; i++) {
		char c = ' ';
	    int query_id = -1;
		fin_query >> c >> query_id >> cnt_node_query >> sum_degree_cur;

		if(i==0)
				stored_emb.resize(num_k);
		if (c != 't'){
			cerr << "No More Query Existed in the Query File!" << endl;
			break;
		}

		mapping_found = 0;

#ifdef WINDOWS
		thread thread_timer(timer, nullptr);
#else
		pthread_t thread_timer;
		pthread_create(&thread_timer, NULL, timer, NULL);
		pthread_detach(thread_timer);
#endif

		TOTAL_BEGIN;

		readQueryGraph();
		coreDecomposition_query();
		isTree = true;
		for (int i = 0; i < cnt_node_query; i ++){
			if (core_number_query[i] >= 2){
				isTree = false;
				break;
			}
        }
        //if(isTree) //TODO: DEBUG wheather it work for tree input
        //    cerr<<"tree"<<endl;
		

        root_node_id = start_node_selection();
        extractResidualStructures();

        BFS_DAG();
			
		topDownInitial();
		bottomUpIterate();
		topDownIterate();
#ifndef MYORDER
		sortByDeg();
#endif
		adjacencyListConstruction();
		//total_candidate_count += getCandidateCount();

        if(use_path_size_order)
            init_weight_by_min();
		// buildGlobalMatchingOrder();
		cur_k = 0;

		SEARCH_BEGIN;
        mapping_found = backtrack_overall();
		SEARCH_END;
		TOTAL_END;

		clearMemory();  //clear indexSet from memory.
#ifdef WINDOWS
		b_search_end = true;
		thread_timer.join();
#else
		pthread_cancel(thread_timer);
#endif

		addUpStatisticParameters(i);

        if(b_over_time_limit){
            cerr<<"TIME LIMIT: "<<i<<"th query"<<endl;
            cerr << "Data File :" << datagraphFile << endl;
            cerr << "Query file:" << querygraphFileFolder << endl;
        }
//        cout << mapping_found << endl;
		cout << i << " " << time_search << " " << time_total << " 0 " << recursive_call_count << " " << mapping_found << " " << overlap_bound << " "<< b_over_time_limit << endl;

		cout << "Coverage "<< computeCoverage() << endl;

		if(b_over_time_limit)
		{
			sum_time_total -= time_total;
			sum_search_time -= time_search;
			sum_mapping -= mapping_found;
		}
        else
		{
			recursive_call_count_sum += recursive_call_count;
			++cnt_end_query;
		}

		b_over_time_limit = false;
#ifdef WINDOWS
		b_search_end = false;
#endif
        recursive_call_count= 0;
	}

	cout << "Average Search Time Per Query => " << setprecision(4) << fixed << sum_search_time / (cnt_end_query) << " ms!" << endl;
	cout << "Average Total Time Per Query => " << setprecision(4) << fixed << sum_time_total / (cnt_end_query) << " ms!" << endl;
	cout << "Total Recursive Call Count => " << recursive_call_count_sum << endl;
	cout << "Total Number of Found Matches=> " << withCommasNoDot(sum_mapping) << endl;
	//cout << "total_candidate_count\t" << total_candidate_count << endl;
	cout << "#End Queries=> " << cnt_end_query << endl;

	fin_query.close();

	return 0;
}
