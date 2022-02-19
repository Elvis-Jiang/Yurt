#pragma once
#include "hashfunction.h"
#include<iostream>
#include<string>
#include<vector>
#include<queue>
#include<set>
#include<map>
#include<cmath>
#include<stdlib.h>
#include<bitset>
#include<memory.h>
#include <sys/time.h>
#include<time.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;

#define UNMAX 0xffffffff
//the threadNum for the base level and the paralleled level num
#define BASETHREADNUM 32
#define PARALLELCOUNTER 2

#define prime 739
#define bigger_p 1048576
#define timer 5
#define M 80000
typedef unsigned int hash_size;
typedef unsigned int weight_size;
typedef unsigned int key_size;
typedef unsigned int matrix_width_size;
typedef unsigned int finger_print_size;
typedef unsigned long long uint64;


struct  hashvalue
{
	unsigned int key;
	int IDnum;
};
bool mycomp(const hashvalue &hv1, const hashvalue &hv2)
{
	return hv1.key<hv2.key;
}
bool operator==(const hashvalue &hv1, const hashvalue &hv2)
{
	return hv1.key==hv2.key;
}
int countjoin( vector<hashvalue> &V1, vector<hashvalue> &V2)
{
 int i1=0,i2=0;
 int count=0;
 while(i1<V1.size())
 {
 	if (i2>=V2.size()) return count;
 	while(V2[i2].key<V1[i1].key)
 	{
 		i2++;
 		if (i2>=V2.size()) return count;
 	}
 	if(V2[i2].key==V1[i1].key)
 	{
 		count+=V1[i1].IDnum;
 	    i1++;
 		i2++;
 		continue;
	 }
	else if(V2[i2].key>V1[i1].key)
	{
 			i1++;
	}
  }
  return count;
}
template<class T>
class hashTableNode
{
public:
	T value;
	unsigned int key;
	hashTableNode<T> *next;
};

template<class T>
class hashTable
{
public:
	hashTableNode<T> **table;
	int tableSize;
	hashTable(int s):tableSize(s)
	{
		table = new hashTableNode<T>*[s];
		memset(table, NULL, tableSize * sizeof(hashTableNode<T>*));
	}
	hashTable()
	{
	}
	void init(int s)
	{
		tableSize = s;
		table = new hashTableNode<T>*[s];
		for(int i=0;i<s;i++)
			table[i]=NULL;
	//	memset(table, NULL, tableSize * sizeof(hashTableNode<T>*));
	}
	~hashTable()
	{
		cleanupHashTable();
		delete [] table;
	}
	void cleanupHashTable()
	{
		hashTableNode<T>*np, *tmp;
		for (int i = 0; i < tableSize; ++i)
		{
			if (table[i] != NULL)
			{
				np = table[i];
				while (np != NULL)
				{
					tmp = np->next;
					delete np;
					np = tmp;
				}
			}
		}
	}
	void insert(unsigned int hash, T value)
	{
		hashTableNode<T> *np;
		bool inTable;
		np = table[hash%tableSize];
		inTable = false;
		for (; np != NULL; np = np->next)
		{
			if (np->key == hash && np->value == value)
			{
				inTable = true;
				break;
			}
		}
		if (!inTable)
		{
			hashTableNode<T>* newNode = new hashTableNode<T>;
			newNode->key = hash;
			newNode->value = value;
			newNode->next = table[hash%tableSize];
			table[hash%tableSize] = newNode;
		}
	}
	void getID(unsigned int hash, vector<T>&IDs)
	{
		hashTableNode<T> *np;
		np=table[hash%tableSize];
		for(;np != NULL; np=np->next)
		{
			if(np->key==hash)
			{
				IDs.push_back(np->value);
			}
		}
		return;
	}
	int countIDnums(unsigned int hash)
{
	int num=0;
	hashTableNode<T> *np;
	np=table[hash%tableSize];
	for(;np!=NULL;np=np->next)
	{
		if(np->key==hash)
			num++;
	}
	return num;
}
};
struct mapnode
{
	unsigned int h;
	unsigned short g;
};

struct linknode
{
	unsigned int key;
	unsigned int weight;
	linknode* next;
};
struct Gbasket
{
	finger_print_size src;
	finger_print_size dst;
	unsigned short idx;
	unsigned int weight;
};
struct LDGbasket{
    unsigned short idx;
	unsigned int weight;
};
struct basketNode{
    LDGbasket *value;
    uint64 *src,*dst;

    vector<linknode*> buffer;
    map<unsigned int, int> index;
    int bufferSize,n;

    basketNode* next[4];
    basketNode(){
        for(int i=0;i<4;i++){
            next[i]=NULL;
        }
        value=NULL;
        bufferSize=0;
        n=0;
    }
};

class GSS
{
private:
    class hashTable<string> mapTable;
	int w;
	int r;
	int p;
	int s;
	int f;
	bool useT;
	int tablesize;
	int bufferSize;

	Gbasket* value;

	public:
		vector<linknode*> buffer;
		map<unsigned int, int> index;
		int n;
		int edge_num; // count the number of edges in the buffer to assist buffer size analysis. Self loop edge is not included as it does not use additional memory.
		GSS(int width, int range, int p_num, int size,int f_num, bool usetable, int tablesize=0);
		~GSS()
		{
			delete[] value;
			cleanupBuffer();

		 }
		 void insert(string s1, string s2,int weight,uint64& insertTime);
		 void cleanupBuffer();
		 int edgeQuery(string s1, string s2);
		 bool query(string s1, string s2);
		 int nodeValueQuery(string s1, int type);//src_type = 0 dst_type = 1
		 int nodeDegreeQuery(string s1, int type);//src_type = 0 dst_type = 1
		 void nodeSuccessorQuery(string s1, vector<string> &IDs);
		 void nodePrecursorQuery(string s2, vector<string> &IDs);
		 int TriangleCounting();
		 void matrixUserate();
};

class LDGSS
{
private:
    class hashTable<string> mapTable;
	int w;
	int r;
	int p;
	int s;
	int f;
	unsigned int max_f_size;
	uint64 matrixSize;
	bool useT;
	int tablesize;
	int maxBufferSize;
	int *finger_size,*finger_per,*get_finger;


	public:
	    basketNode *valueTree;
		int edge_num; // count the number of edges in the buffer to assist buffer size analysis. Self loop edge is not included as it does not use additional memory.
		LDGSS(int width, int range, int p_num, int size,int f_num, double bufferRate, bool usetable, int tablesize=0);
		~LDGSS()
		{
			delete[] valueTree;
			cleanupBuffer();

		 }
		 basketNode* getBasketNode(unsigned int fp1,unsigned int fp2,int& level,uint64& pre);
		 void insert(string s1, string s2,int weight,uint64& insertTime);
		 void cleanupBuffer();
		 int edgeQuery(string s1, string s2);
		 bool query(string s1, string s2);
		 int nodeValueQuery(string s1, int type);//src_type = 0 dst_type = 1
		 int nodeDegreeQuery(string s1, int type);//src_type = 0 dst_type = 1
		 void nodeSuccessorQuery(string s1, vector<string> &IDs);
		 void nodePrecursorQuery(string s2, vector<string> &IDs);
		 int TriangleCounting();
		 double matrixUserate(basketNode *basketCurr);
		 bool appendBasketTree(basketNode* basketCurr,int level);
		 void getSuccessorInTree(basketNode * basketCurr,int *tmp1,matrix_width_size h1,finger_print_size preg1,finger_print_size g1,
                        finger_print_size pre,finger_print_size level,vector<string>&IDs,int type,unsigned int& res);
		 void getPrecursorInTree(basketNode * basketCurr,int *tmp1,matrix_width_size h1,finger_print_size preg1,finger_print_size g1,
                           finger_print_size pre,finger_print_size level,vector<string>&IDs,int type,unsigned int& res);
		 void insertBasket(weight_size weight,key_size k1,key_size k2,finger_print_size g1,finger_print_size g2,
                     basketNode* basketCurr,int finger_sizeCurr,int finger_perCur,bool append,int level=0);
		 basketNode* initialBasketNode(int need64);
};
int calHash(int tmp_g,int tmp_s,int k,int w){
    int shifter = tmp_g;
    for (int v = 0; v < tmp_s; v++)
        shifter = (shifter*timer + prime) % bigger_p;
    int tmp_h = k;
    while (tmp_h < shifter)
        tmp_h += w;
    tmp_h -= shifter;
    return tmp_h;
}

struct edgeval{
    unsigned int val;
    string from,to;
    edgeval(string _from,string _to,unsigned int _val){
        from=_from;to=_to;val=_val;
    }
};



bool operator < (const edgeval& a,const edgeval& b)
{
    return a.val<b.val;
}

class heap{
private:
    //将以k为根节点的子树设置成小根堆
    //其中左右子树均已为小根堆
    void heapSubTree(int k,int n){
        while((2*k+1)<=n){
            int temp=2*k+1;
            if((2*k+2)<=n&&(*num[2*k+2])<(*num[2*k+1])){
                temp=2*k+2;
            }
            if((*num[temp])<(*num[k])){
                edgeval *t=num[temp];
                num[temp]=num[k];
                num[k]=t;
                k=temp;
            }
            else break;
        }
    }
public:
    edgeval **num;
    int s;
    heap(){}
    //将小根堆初始化
    void initialHeap(int l){
        s=l;
        num=new edgeval*[s];
        for(int i=0;i<s;i++)num[i]=new edgeval("null","null",0);
        for(int j=(s-2)/2;j>=0;j--){
            heapSubTree(j,s-1);
        }
    }
//    bool findincrease(string s1,string s2,int h,unsigned int edgevalue){
//        unsigned int h1=(*hfunc[h])((unsigned char*)(s1.c_str()), s1.length());
//        unsigned int h2=(*hfunc[h])((unsigned char*)(s2.c_str()), s2.length());
//        for(int i=0;i<s;i++){
//            unsigned int h11=(*hfunc[h])((unsigned char*)(num[i]->from.c_str()), num[i]->from.length());
//            unsigned int h22=(*hfunc[h])((unsigned char*)(num[i]->to.c_str()), num[i]->to.length());
//            if(h1==h11&&h2==h22){
//                if(edgevalue-num[i]->val>num[0]->val){
//                    return false;
//                }
//                else{
//                    num[i]->from=s1;
//                    num[i]->to=s2;
//                    num[i]->val=edgevalue;
//                    heapSubTree(i,s-1);
//                    cout<<"find and increase\n";
//                    return true;
//                }
//            }
//        }
//        return false;
//    }
    unsigned int find(string s1,string s2){
        for(int i=0;i<s;i++){
            if(num[i]->from==s1&&num[i]->to==s2){
                return num[i]->val;
            }
        }
        return -1;
    }
    ~heap(){
        for(int i=0;i<s;i++){
            delete num[i];
        }
    }
    void pop(){
        delete num[0];
        num[0]=num[s-1];
        s--;
        heapSubTree(0,s-1);
    }
    void push(string s1,string s2,unsigned int weight){
        num[s]=new edgeval(s1,s2,weight);s++;
        int temp=(s-2)/2,child=s-1;
        while(temp>=0){
            if((*num[child])<(*num[temp])){
                edgeval *t=num[child];
                num[child]=num[temp];
                num[temp]=t;
                child=temp;
                temp=(temp-1)/2;
            }
            else break;
        }
    }
    bool increase(string s1,string s2,unsigned int weight){
        for(int i=0;i<s;i++){
            if(num[i]->from==s1&&num[i]->to==s2){
                int now=i;
                num[now]->val=weight;
                heapSubTree(now,s-1);
                return true;
            }
        }
        return false;
    }
    unsigned int top(){
        return num[0]->val;
    }
    void topres(string &s1,string &s2){
        s1=num[0]->from;
        s2=num[0]->to;
    }
    bool empty(){
        if(!s)return true;
        return false;
    }
    int size(){
        return s;
    }

};
