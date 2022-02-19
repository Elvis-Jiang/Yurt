#ifndef GSS_H_INCLUDED
#define GSS_H_INCLUDED
#include "querysupportstruct.h"
using namespace std;

GSS::GSS(int width, int range, int p_num, int size,int f_num, bool usehashtable, int TableSize)
//the side length of matrix, the length of hash address list, the number of candidate bucket
// the number of rooms, whether to use hash table
//and the size of the table.
// Hash table which stores the original nodes can be omitted if not needed. For node query,
//  reachability, edge query not needed. But needed for triangle counting, degree query, and successor / precursor queries.
{
    bufferSize=0;
	w = width;
	r = range; /* r x r mapped baskets */
	p = p_num; /*candidate buckets*/
	s = size; /*multiple rooms*/
	f = f_num; /*finger print length*/
	n = 0;
	value = new Gbasket[w*w];
	useT=usehashtable;
	tablesize=TableSize;
	memset(value, 0, sizeof(Gbasket)*w*w);
	if(usehashtable)
	mapTable.init(tablesize);
}
void GSS::matrixUserate(){
    double used=0;
    for(int i=0;i<w*w;i++){
        {
            if(value[i].src)used+=1;
        }
    }
    cout<<bufferSize<<" "<<used/(w*w)<<endl;
}
void GSS::cleanupBuffer()
{
	vector<linknode*>::iterator IT = buffer.begin();
	linknode* e, *tmp;
	for (; IT != buffer.end(); ++IT)
	{
		e = *IT;
		while (e != NULL)
		{
			tmp = e->next;
			delete e;
			e = tmp;
		}
	}
}
void GSS::insert(string s1, string s2,int weight,uint64& insertTime)// s1 is the ID of the source node, s2 is the ID of the destination node, weight is the edge weight.
{
		timeval t_start, t_end;
        gettimeofday( &t_start, NULL);
		unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
		unsigned int hash2 = (*hfunc[0])((unsigned char*)(s2.c_str()), s2.length());
		unsigned int tmp = pow(2,f)-1;
		//取低f位
		unsigned short g1 = hash1 & tmp;
		if(g1==0) g1+=1;
		unsigned int h1 = (hash1>>f)%w;
		unsigned short g2 = hash2 & tmp;
		if(g2==0) g2+=1;
		unsigned int h2 = (hash2>>f)%w;
		//到时候还原到k
		unsigned int k1 = (h1<<f)+g1;
		unsigned int k2 = (h2<<f)+g2;

		if(useT){
		mapTable.insert(k1, s1);
		mapTable.insert(k2, s2);
		}

		int* tmp1 = new int[r];
		int* tmp2 = new int[r];
		tmp1[0] = g1;
		tmp2[0] = g2;
		for(int i=1;i<r;i++)
		{
			tmp1[i]=(tmp1[i-1]*timer+prime)%bigger_p;
			tmp2[i]=(tmp2[i-1]*timer+prime)%bigger_p;
		}
		bool inserted=false;
		long key = g1+g2;
		for(int i=0;i<p;i++)
		{
			key = (key*timer+prime)%bigger_p;
			int index = key%(r*r);
			int index1 = index/r;
			int index2 = index%r;
			int p1 = (h1+tmp1[index1])%w;
			int p2 = (h2+tmp2[index2])%w;


			int pos = p1*w + p2;
			{
				if ( ( value[pos].idx==(index2|(index1<<8)) ) && (value[pos].src== g1) && (value[pos].dst == g2) )
				{
					value[pos].weight += weight;
					inserted = true;
					break;
				}
				if (value[pos].src == 0)
				{
					value[pos].idx = (index2|(index1<<8));
					value[pos].src = g1;
					value[pos].dst = g2;
					value[pos].weight = weight;
					inserted = true;
					break;
				}
			}
			if(inserted)
				break;
		}
		if(!inserted)
		{

			map<unsigned int, int>::iterator it = index.find(k1);
			if(it!=index.end())
			{
				int tag = it->second;
				linknode* node = buffer[tag];
				while(true)
				{
					if (node->key == k2)
					{
						node->weight += weight;
						break;
					}
					if(node->next==NULL)
					{
						linknode* ins = new linknode;
						bufferSize++;
						ins->key = k2;
						ins->weight = weight;
						ins->next = NULL;
						node->next = ins;
						break;
					}
					node = node->next;
				}
			}
			else
			{
				index[k1] = n;
				n++;
				linknode* node = new linknode;
				bufferSize++;
				node->key = k1;
				node->weight = 0;
				if (k1 != k2)//k1==k2 means loop
				{
					linknode* ins = new linknode;
					bufferSize++;
					ins->key = k2;
					ins->weight = weight;
					ins->next = NULL;
					node->next = ins;
				}
				else
				{
					node->weight += weight;
					node->next = NULL;
				}
				buffer.push_back(node);
			}
		}
		delete [] tmp1;
		delete [] tmp2;
		gettimeofday( &t_end, NULL);
        insertTime+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;
	return;
}
void GSS::nodeSuccessorQuery(string s1, vector<string>&IDs)// query the successors of a node, s1 is the ID of the queried node. results are put in the vector, hash table needed.
{
	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
	int tmp=pow(2,f)-1;
	unsigned short g1=hash1 & tmp;
	if(g1==0) g1+=1;
	unsigned int h1 = (hash1>>f)%w;
	unsigned int k1 = (h1 << f) + g1;
	int* tmp1 = new int[r];
	tmp1[0] = g1;
	for (int i = 1; i < r; i++)
	{
		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
	}
	for (int i = 0; i < r; i++)
	{
		int p1 = (h1 + tmp1[i]) % w;
		for (int k = 0; k < w; k++)
		{
			int pos = p1*w + k;
			{
				if ((value[pos].idx>>8) == i && (value[pos].src == g1))
				{
					     int tmp_g = value[pos].dst;
						 int tmp_s = (value[pos].idx&((1<<8)-1));
						 int val=(calHash(tmp_g,tmp_s,k,w)<<f)+tmp_g;
						 mapTable.getID(val, IDs);
				}
			}
		}
	}
		map<unsigned int, int>::iterator it = index.find(k1);
		if (it != index.end())
		{
			int tag = it->second;
			linknode* node = buffer[tag];
			//node = node->next;
			while (node != NULL)
			{
			    if(node->weight)
                    mapTable.getID(node->key, IDs);
				node=node->next;
			}
		}
		delete []tmp1;
		return;
}
void GSS::nodePrecursorQuery(string s1, vector<string>&IDs) // same as successor query
{
	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
	int tmp=pow(2,f)-1;
	unsigned short g1=hash1 & tmp;
	unsigned int h1 = (hash1>>f)%w;
	if(g1==0) g1+=1;
	int* tmp1 = new int[r];
	tmp1[0] = g1;
	unsigned int k1 = (h1 << f) + g1;
	for (int i = 1; i < r; i++)
	{
		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
	}
	for (int i = 0; i < r; i++)
	{
		int p1 = (h1 + tmp1[i]) % w;
		for (int k = 0; k < w; k++)
		{
			int pos = p1 + k*w;
			{
				if ((value[pos].idx&((1<<8)-1)) == i && (value[pos].dst == g1))
				{
					     int tmp_g = value[pos].src;
						 int tmp_s = (value[pos].idx>>8);

						 int val=(calHash(tmp_g,tmp_s,k,w)<<f)+tmp_g;
						 mapTable.getID(val, IDs);
				}
			}
		}
	}
			for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
		{
			int tag = it->second;
			int src = it->first;
			linknode* node = buffer[tag];
			//node = node->next;
			while (node != NULL)
			{
				if(node->key == k1&&node->weight)
				{
					mapTable.getID(src, IDs);
					break;
				}
				node = node->next;
			}
		}
		delete []tmp1;
		return;
}
int GSS::edgeQuery(string s1, string s2)// s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
{
	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
	unsigned int hash2 = (*hfunc[0])((unsigned char*)(s2.c_str()), s2.length());
	int tmp = pow(2, f) - 1;
	unsigned short g1 = hash1 & tmp;
	if (g1 == 0) g1 += 1;
	unsigned int h1 = (hash1 >> f) % w;
	unsigned short g2 = hash2 & tmp;
	if (g2 == 0) g2 += 1;
	unsigned int h2 = (hash2 >> f) % w;
	int* tmp1 = new int[r];
	int* tmp2 = new int[r];
	tmp1[0] = g1;
	tmp2[0] = g2;
	for (int i = 1; i<r; i++)
	{
		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
		tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
	}
	long key = g1 + g2;

	for (int i = 0; i<p; i++)
	{
		key = (key * timer + prime) % bigger_p;
		int index = key % (r*r);
		int index1 = index / r;
		int index2 = index%r;
		int p1 = (h1 + tmp1[index1]) % w;
		int p2 = (h2 + tmp2[index2]) % w;
		int pos = p1*w + p2;
		{

			if ((value[pos].idx == (index2 | (index1 << 8))) && (value[pos].src == g1) && (value[pos].dst == g2))
			{
				delete []tmp1;
				delete []tmp2;
				return value[pos].weight;
			}
		}

	}
	unsigned int k1 = (h1 << f) + g1;
	unsigned int k2 = (h2 << f) + g2;
	map<unsigned int, int>::iterator it = index.find(k1);
	if (it != index.end())
	{
		int tag = it->second;
		linknode* node = buffer[tag];
		while (node!=NULL)
		{
			if (node->key == k2)
			{
				delete []tmp1;
				delete []tmp2;
				return node->weight;
			}
			node = node->next;
		}
	}
	//cout<<"could not find the edge\n";
    delete []tmp1;
    delete []tmp2;
    return 0;
}
/*type 0 is for successor query, type 1 is for preccessor query*/
int GSS::nodeValueQuery(string s1, int type) // s1 is the ID of the queried node, function for node query.
{
	int weight = 0;
	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
	int tmp = pow(2, f) - 1;
	unsigned short g1 = hash1 & tmp;
	if (g1 == 0) g1 += 1;
	unsigned int h1 = (hash1 >> f) % w;
	int* tmp1 = new int[r];
	tmp1[0] = g1;
	for (int i = 1; i < r; i++)
	{
		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
	}
	for (int i = 0; i < r; i++)
	{
		int p1 = (h1 + tmp1[i]) % w;
		for (int k = 0; k < w; k++)
		{
			if (type == 0)/*successor query*/
			{
				int pos = p1*w + k;
				{
					if ((value[pos].idx>>8) == i && (value[pos].src == g1))
					{
						weight += value[pos].weight;
					}
				}
			}
			else if (type == 1)/*precursor query*/
			{
				int pos = p1 + k*w;
				{
					if ((value[pos].idx&((1<<8)-1)) == i && (value[pos].dst == g1))
					{
						weight += value[pos].weight;
					}
				}
			}
		}
	}
	if (type == 0)
	{
		unsigned int k1 = (h1 << f) + g1;
		map<unsigned int, int>::iterator it = index.find(k1);
		if (it != index.end())
		{
			int tag = it->second;
			linknode* node = buffer[tag];
			//node = node->next;
			while (node != NULL)
			{
				weight += node->weight;
				node = node->next;
			}
		}
	}
	else if (type==1)
	{
		unsigned int k1 = (h1 << f) + g1;
		for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
		{
			int tag = it->second;
			linknode* node = buffer[tag];
			//node = node->next;
			while (node != NULL)
			{
				if(node->key == k1)
					weight += node->weight;
				node = node->next;
			}
		}
	}
	delete []tmp1;
	return weight;
}
/*type 0 is for successor query, type 1 is for precusor query*/
int GSS::nodeDegreeQuery(string s1, int type) // s1 is the ID of the queried node, return the in/out degree
{
	int degree = 0;
	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
	int tmp = pow(2, f) - 1;
	unsigned short g1 = hash1 & tmp;
	if (g1 == 0) g1 += 1;
	unsigned int h1 = (hash1 >> f) % w;
	int* tmp1 = new int[r];
	tmp1[0] = g1;
	for (int i = 1; i < r; i++)
	{
		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
	}
	for (int i = 0; i < r; i++)
	{
		int p1 = (h1 + tmp1[i]) % w;
		for (int k = 0; k < w; k++)
		{
			if (type == 0)/*successor query*/
			{
				int pos = p1*w + k;
				{
					if ((value[pos].idx>>8) == i && (value[pos].src == g1))
					{
						 int tmp_g = value[pos].dst;
						 int tmp_s = (value[pos].idx&((1<<8)-1));

						 int val=(calHash(tmp_g,tmp_s,k,w)<<f)+tmp_g;


						 degree+=mapTable.countIDnums(val);
					}
				}
			}
			else if (type == 1)/*precursor query*/
			{
				int pos = p1 + k*w;
				{
					if ((value[pos].idx&((1<<8)-1)) == i && (value[pos].dst == g1))
					{
						 int tmp_g = value[pos].src;
						 int tmp_s = (value[pos].idx>>8);

						 int val=(calHash(tmp_g,tmp_s,k,w)<<f)+tmp_g;

						degree+=mapTable.countIDnums(val);
					}
				}
			}
		}
	}

	if (type == 0)
	{
		unsigned int k1 = (h1 << f) + g1;
		map<unsigned int, int>::iterator it = index.find(k1);
		if (it != index.end())
		{
			int tag = it->second;
			linknode* node = buffer[tag];
			//node = node->next;
			while (node != NULL)
			{
			    if(node->weight)
                    degree+=mapTable.countIDnums(node->key);
				node = node->next;
			}
		}
	}
	else if (type == 1)
	{
		unsigned int k1 = (h1 << f) + g1;
		for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
		{
			int tag = it->second;
			linknode* node = buffer[tag];
			unsigned int src=node->key;
			//node = node->next;
			while (node != NULL)
			{
				if (node->key == k1)
				{
				    if(node->weight)
                        degree+=mapTable.countIDnums(src);
					break;
				 }
				node = node->next;
			}
		}
	}
	delete[]tmp1;
	return degree;
}
#endif // GSS_H_INCLUDED
