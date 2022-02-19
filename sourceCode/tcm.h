#ifndef TCM_H
#define TCM_H
#include<math.h>
#include<string>
#include<iostream>
#include<memory.h>
#include<queue>
#include <fstream>
#include <bitset>
#include <set>
#include "querysupportstruct.h"
#include "hashTable.h"
#include <sys/time.h>
typedef unsigned int AL;
using namespace std;

class tcm
{
private:
    int **widths;
    unsigned int hashnum;
    AL **value;
    int w;

    BOBHash64 * bobhash[MAX_HASH_NUM];
public:
    tcm(unsigned w,unsigned int h_num);
    ~tcm()
    {
        for(int i=0;i<hashnum;i++){
            delete [] widths[i];
            delete [] value[i];
        }
        delete [] value;
        delete [] widths;
        delete [] bobhash;
    }
    void insert(string s1, string s2,unsigned int weight, unsigned int& hashTime, unsigned int& updateTime);
    unsigned int edgeQuery(string s1, string s2);
    unsigned int nodeQuery(string s1, bool type);//outval=1,inval=0
    bool reachbility(string s1,string s2);
    void print(int num,int cou);
    void bitUseRate();
    void useRate();

};
//void tcm::print(int num,int cou){
//    ofstream location_out;
//    string str="./dataDis/lkml2/"+to_string(cou)+".txt";
//    char*ps=(char*)str.data();
//    location_out.open(ps, std::ios::out | std::ios::app);
//    for(int i=0;i<widths[num][0];i++){
//        for(int j=0;j<widths[num][1];j++){
//            location_out<<value[0][i*widths[num][1]+j]<<"\n";
//        }
//    }
//
//}

void tcm::useRate(){
    double useRate=0.0;
    for(int i=0;i<hashnum;i++){
        double mid=0.0;
        int cou=widths[i][0]*widths[i][1];
        for(int j=0;j<cou;j++){
            if(value[i][j])mid+=1;
        }
        useRate+=mid/cou;
    }
    cout<<useRate/hashnum<<endl;
}

void tcm::bitUseRate(){
    double useRate=0.0;
    for(int i=0;i<hashnum;i++){
        double mid=0.0;
        int cou=widths[i][0]*widths[i][1];
        for(int j=0;j<cou;j++){
            if(value[i][j])mid+=(log(value[i][j]+1)/log(2))/16;
        }
        useRate+=mid/cou;
    }
    cout<<useRate/hashnum<<endl;
}
unsigned int tcm::nodeQuery(string s1, bool type){
    unsigned int m=UNMAX;
	for(int i = 0; i < hashnum; i++){
        unsigned int inner=0;
		if(type){
            if(i%2)continue;
            unsigned int hash1;
            hash1=(bobhash[i]->run(s1.c_str(), strlen(s1.c_str())))%widths[i][0];
            for(int j=0;j<widths[i][1];j++){
                inner+=value[i][hash1*widths[i][1]+j];
            }
		}
		else{
		    if(i!=0&&i%2==0)continue;
		    unsigned int hash1;
            hash1=(bobhash[i]->run(s1.c_str(), strlen(s1.c_str())))%widths[i][1];
            for(int j=0;j<widths[i][0];j++){
                inner+=value[i][hash1+j*widths[i][1]];
            }
		}
		if(inner<m){
            m=inner;
            if(m==0)return 0;
		}
    }
    if(m==UNMAX)return 0;
    return m;
}
tcm::tcm(unsigned w,unsigned int h_num){
    hashnum=h_num;
    widths=new int*[h_num];
    for(int i=0;i<h_num;i++)widths[i]=new int[2];
    widths[0][0]=sqrt(w);
    this->w=w;
    widths[0][1]=widths[0][0];
    for(int i=1;i<=(h_num>>1);i++){
        int divd=1<<(2*i);
        widths[2*i-1][0]=sqrt(w)/divd;
        if(widths[2*i-1][0]==0)widths[2*i-1][0]=1;
        widths[2*i-1][1]=w/widths[2*i-1][0];
//        widths[2*i-1][1]+=(widths[2*i-1][0]-(widths[2*i-1][1]%widths[2*i-1][0]));

        widths[2*i][0]=widths[2*i-1][1];
        widths[2*i][1]=widths[2*i-1][0];
//        widths[2*i][0]+=(widths[2*i][1]-(widths[2*i][0]%widths[2*i][1]));
    }
    value=new AL*[hashnum];
    for(int i=0;i<hashnum;i++){
        int counternum=widths[i][0]*widths[i][1];
        value[i]=new AL[counternum];
        memset(value[i],0,sizeof(AL)*counternum);
    }

    for(int i = 0; i < MAX_HASH_NUM; i++){
		bobhash[i] = new BOBHash64(i + 1000);
	}
}

//广度优先搜索探寻两点之间是否可达
bool tcm::reachbility(string s1,string s2){
    if(s1==s2)return true;
    for(int i=0;i<hashnum;i++){
        set<unsigned int> visited;
        queue<unsigned int> nei;
        bool reach=false;
        unsigned int hash1= (bobhash[i]->run(s1.c_str(), strlen(s1.c_str())));
		unsigned int hash2 = (bobhash[i]->run(s2.c_str(), strlen(s2.c_str())));
		if(widths[i][0]<widths[i][1]){
            unsigned int h1=hash1%widths[i][1],h2=hash2%widths[i][1];
            visited.insert(h1);
            nei.push(h1);
            while(!nei.empty()&&!reach){
                unsigned int now=nei.front()%widths[i][0];
                nei.pop();
                for(int j=0;j<widths[i][1];j++){
                    if(value[i][now*widths[i][1]+j]){
                        if(j==h2){
                            reach=true;
                            break;
                        }
                        else if(!visited.count(j)){
                            visited.insert(j);
                            nei.push(j);
                        }
                    }

                }
            }
		}
		else{
		    unsigned int mag=widths[i][0]/widths[i][1];
            unsigned int h1=hash1%widths[i][0],h2=hash2%widths[i][1];
            visited.insert(h1);
            nei.push(h1);
            while(!nei.empty()&&!reach){
                unsigned int now=nei.front();
                nei.pop();
                for(int j=0;j<widths[i][1];j++){
                    if(value[i][now*widths[i][1]+j]){
                        if(j==h2){
                            reach=true;
                            break;
                        }
                        else if(!visited.count(j)){
                            for(int k=0;k<mag;k++){
                                visited.insert(k*widths[i][1]+j);
                                nei.push(k*widths[i][1]+j);
                            }
                        }
                    }

                }
            }
		}
        if(!reach)return false;
    }
    return true;
}

unsigned int tcm::edgeQuery(string s1, string s2){
    unsigned int m=UNMAX;
    unsigned hash1,hash2;
	for (int i = 0; i < hashnum; i++){
        hash1=(bobhash[i]->run(s1.c_str(), strlen(s1.c_str())))%widths[i][0];
        hash2=(bobhash[i]->run(s2.c_str(), strlen(s2.c_str())))%widths[i][1];
		if(value[i][hash1*widths[i][1]+hash2]<m){
            m=value[i][hash1*widths[i][1]+hash2];
            if(m==0)return m;
		}
    }
    if(m==UNMAX)return 0;
    return m;
}

void tcm::insert(string s1, string s2,unsigned int weight, unsigned int& hashTime, unsigned int& updateTime){
    timeval t_start, t_end;
    gettimeofday( &t_start, NULL);
    unsigned hash1,hash2;
    for(int i=0;i<hashnum;i++){
        hash1=(bobhash[i]->run(s1.c_str(), strlen(s1.c_str())))%widths[i][0];
        hash2=(bobhash[i]->run(s2.c_str(), strlen(s2.c_str())))%widths[i][1];
		//cout<<"a";
		value[i][hash1*widths[i][1]+hash2]+=weight;
    }
    gettimeofday( &t_end, NULL);
    hashTime+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;
}
#endif // TCM_H
