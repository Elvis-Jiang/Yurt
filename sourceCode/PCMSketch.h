#ifndef _PCMSKETCH_H
#define _PCMSKETCH_H

#include <algorithm>
#include <cstring>
#include <string.h>
#include "querysupportstruct.h"
#include <iostream>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
using namespace std;

typedef unsigned long long int uint64;


void* PCMSketch::pAppend(void* args){
    threadParameter *tp;
    tp = (threadParameter*)(args);
    int num=tp->num,level=tp->level,from=tp->from,to=tp->to,type=tp->type;
    PCMSketch *pcm=tp->pthis;
    if(type==0){
        for(int pos=from;pos<to;pos++){
            int my_word_index=(pos/pcm->per_base),counter_offset=(pos%pcm->per_base);
            uint64 downValue=((pcm->allcounter[num][level][my_word_index] >> (counter_offset*pcm->base_counter_size)) & pcm->max_base);

            int word_index=(pos/64),offset=(pos%64);
            uint64 left=(downValue>>1),value=(downValue&0x1);

            pcm->allcounter[num][level-1][word_index]|= (value << offset);
            pcm->accessedFlag[num][level-pcm->base_level[num]][pos]=true;
        }
    }

    else{
        int per,counter_size,max_here;
        if(level==pcm->base_level[num]+1){
            per=pcm->per_base;
            counter_size=pcm->base_counter_size;
        }
        else{
            per=pcm->per_magnify;
            counter_size=pcm->magnify_counter_size;
        }
        max_here=(1<<(counter_size-1))-1;

        for(int pos=from;pos<to;pos++){
            int my_word_index=(pos/pcm->per_magnify),counter_offset=(pos%pcm->per_magnify);
            uint64 downValue=((pcm->allcounter[num][level][my_word_index] >> (counter_offset*pcm->magnify_counter_size)) & pcm->max_magnify);
            bool overflowFlag=((pcm->allcounter[num][level][my_word_index] >>
                                (counter_offset*pcm->magnify_counter_size)) & (uint64)(1<<(pcm->magnify_counter_size-1)));

            uint64 left=(downValue>>1),value=(downValue&0x1);


            bool overflowed=(overflowFlag||left);

            vector<int> downPos;
            pcm->get_down_pos(pos,level,num,downPos);


            for(int z=0;z<downPos.size();z++){
                int b=downPos[z];
                int word_indexdown=(b/per),offsetdown=(b%per);
                while(!pcm->accessedFlag[num][level-pcm->base_level[num]-1][b]){}

                uint64 downValue=((pcm->allcounter[num][level-1][word_indexdown] >>
                                   (offsetdown*counter_size)) & (max_here));
                downValue>>=1;
                uint64 downValue1=downValue+(value<<(counter_size-2));

                if(((pcm->allcounter[num][level-1][word_indexdown] >>
                                (offsetdown*counter_size)) & (uint64)(max_here+1))){

                    pcm->allcounter[num][level-1][word_indexdown]&=(~((uint64)(max_here) << (offsetdown*counter_size)));
                    pcm->allcounter[num][level-1][word_indexdown]|= (downValue1<<(offsetdown*counter_size));
                    if(!overflowed){
                        pcm->allcounter[num][level-1][word_indexdown]&= (~((uint64)(max_here+1) << (offsetdown*counter_size)));
                        if(level==pcm->base_level[num]+1)pcm->overflow_count[num]-=1;
                    }
                }
                else{
                    pcm->allcounter[num][level-1][word_indexdown]&=(~((uint64)(max_here) << (offsetdown*counter_size)));
                    pcm->allcounter[num][level-1][word_indexdown]|= (downValue<<(offsetdown*counter_size));
                }


                pcm->accessedFlag[num][level-pcm->base_level[num]-1][b]=false;
            }
            if(type==1){
               pcm->accessedFlag[num][level-pcm->base_level[num]][pos]=true;
            }
//            else{
//               pcm->allcounter[num][level][my_word_index]&=(~((uint64)pcm->max_magnify << (counter_offset*pcm->magnify_counter_size)));
//               pcm->allcounter[num][level][my_word_index]|= (left<<(counter_offset*pcm->magnify_counter_size));
//            }
        }
    }

    return NULL;
}
//append the structure
void PCMSketch::down_flow(int num){
    timeval t_start, t_end;
    gettimeofday( &t_start, NULL);
    appendTime+=1;
    uint64 *a=new uint64[wordNum];

    memset(a,0,sizeof(uint64)*wordNum);

    allcounter[num].insert(allcounter[num].begin()+base_level[num],a);
    int dep=depth[num][0],wid=width[num][0];
    depth[num].insert(depth[num].begin()+base_level[num],dep);
    width[num].insert(width[num].begin()+base_level[num],wid);

    base_level[num]++;

    //parallel down flow
    pthread_t *tids[PARALLELCOUNTER];

    threadParameter *tp[PARALLELCOUNTER];

    int threadNum=BASETHREADNUM;
    int type;

    for(int i=base_level[num];i<base_level[num]+PARALLELCOUNTER;i++){//down flow
        tids[i]=new pthread_t[threadNum];
        tp[i]=new threadParameter[threadNum];
        if(i==base_level[num])type=0;
        else if(i==base_level[num]+PARALLELCOUNTER-1)type=2;
        else type=1;
        int counter_num=width[num][i]*depth[num][i],derta=counter_num/threadNum;

        if(derta==0){
            tp[i][0]=threadParameter(this,num,i,0,counter_num,type);

            int ret = pthread_create(&tids[i][0], NULL, pAppend, (void *)&(tp[i][0]));
        }
        else{
            int from=0;
            for(int j=0;j<threadNum;j++){
                int to;
                if(j+1==threadNum)to=counter_num;
                else to=from+derta;
                tp[i][j]=threadParameter(this,num,i,from,to,type);

                int ret = pthread_create(&tids[i][j], NULL, pAppend, (void *)&(tp[i][j]));

                from=to;
            }
        }
        if(threadNum>=4)threadNum>>=2;
        else threadNum=1;
    }

    threadNum=BASETHREADNUM;
    for(int i=base_level[num];i<base_level[num]+PARALLELCOUNTER;i++){//down flow
        int counter_num=width[num][i]*depth[num][i],derta=counter_num/threadNum;
        if(derta==0){
            pthread_join(tids[i][0],NULL);
        }
        else{
            for(int j=0;j<threadNum;j++){
                pthread_join(tids[i][j],NULL);
            }
        }
        if(threadNum>=4)threadNum>>=2;
        else threadNum=1;
    }

    for(int i=base_level[num]+PARALLELCOUNTER;i<allcounter[num].size();i++){
        int level=i;
        if(i==base_level[num])type=0;
        else if(i==allcounter[num].size()-1)type=2;
        else type=1;
        int counter_num=width[num][i]*depth[num][i];
        if(i==base_level[num]){
            for(int pos=0;pos<counter_num;pos++){
                int my_word_index=(pos/per_base),counter_offset=(pos%per_base);
                uint64 downValue=((allcounter[num][level][my_word_index] >> (counter_offset*base_counter_size)) & max_base);

                int word_index=(pos/64),offset=(pos%64);
                uint64 left=(downValue>>1),value=(downValue&0x1);

                allcounter[num][level-1][word_index]|= (value << offset);
            }
        }

        else{
            int per,counter_size,max_here;
            if(level==base_level[num]+1){
                per=per_base;
                counter_size=base_counter_size;
            }
            else{
                per=per_magnify;
                counter_size=magnify_counter_size;
            }
            max_here=(1<<(counter_size-1))-1;

            for(int pos=0;pos<counter_num;pos++){
                int my_word_index=(pos/per_magnify),counter_offset=(pos%per_magnify);
                uint64 downValue=((allcounter[num][level][my_word_index] >> (counter_offset*magnify_counter_size)) & max_magnify);
                bool overflowFlag=((allcounter[num][level][my_word_index] >>
                        (counter_offset*magnify_counter_size)) & (uint64)(1<<(magnify_counter_size-1)));

                uint64 left=(downValue>>1),value=(downValue&0x1);


                bool overflowed=(overflowFlag||left);

                vector<int> downPos;
                get_down_pos(pos,level,num,downPos);


                for(int z=0;z<downPos.size();z++){
                    int b=downPos[z];
                    int word_indexdown=(b/per),offsetdown=(b%per);

                    uint64 downValue=((allcounter[num][level-1][word_indexdown] >>
                            (offsetdown*counter_size)) & (max_here));
                    downValue>>=1;
                    uint64 downValue1=downValue+(value<<(counter_size-2));

                    if(((allcounter[num][level-1][word_indexdown] >>
                            (offsetdown*counter_size)) & (uint64)(max_here+1))){

                        allcounter[num][level-1][word_indexdown]&=(~((uint64)(max_here) << (offsetdown*counter_size)));
                        allcounter[num][level-1][word_indexdown]|= (downValue1<<(offsetdown*counter_size));
                        if(!overflowed){
                            allcounter[num][level-1][word_indexdown]&= (~((uint64)(max_here+1) << (offsetdown*counter_size)));
                            if(level==base_level[num]+1)overflow_count[num]-=1;
                        }
                    }
                    else{
                        allcounter[num][level-1][word_indexdown]&=(~((uint64)(max_here) << (offsetdown*counter_size)));
                        allcounter[num][level-1][word_indexdown]|= (downValue<<(offsetdown*counter_size));
                    }
                }
                if(type==2){
                   allcounter[num][level][my_word_index]&=(~((uint64)max_magnify << (counter_offset*magnify_counter_size)));
                   allcounter[num][level][my_word_index]|= (left<<(counter_offset*magnify_counter_size));
                }
            }
        }
    }


    leftRatio+=overflow_count[num]/(double)up_bound;
    if((int)appendTime%(int)d==0){
        cout<<leftRatio/d<<endl;
        leftRatio=0;
    }
    gettimeofday( &t_end, NULL);
    appendDuration+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;
}
//calculate the overflow rate of the first layer
void PCMSketch::overflowRate(){
//    ofstream location_out;
//    string str=to_string(w)+"overflow.txt";
//    char*ps=(char*)str.data();
//    location_out.open(ps, std::ios::out | std::ios::app);
//    int used=0;
//    for(int k=0;k<depth[0][0]*width[0][0];k++){
//        int word_index,offset,temp;
//        word_index=(k/per_base),offset=k%per_base;
//        temp=base_counter_size;
//
//        if((allcounter[0][0][word_index] >> (offset*temp)) & (uint64)(1<<(temp-1))){
//            used++;
//        }
//        if((k+1)%width[0][0]==0){
//            location_out<<used<<"\n";
//            used=0;
//        }
//    }
//    location_out.close();

    int cou[15];
    double rate[15];
    for(int i=0;i<15;i++){
        cou[i]=0;
        rate[i]=0.0;
    }
    int dep,wid,used;

    for(int i=0;i<d;i++){
        for(int j=0;j<allcounter[i].size();j++){
            dep=depth[i][j];wid=width[i][j];
            used=0;
            for(int k=0;k<dep*wid;k++){
                int word_index,offset,temp;
                if(j==0){
                    word_index=(k/per_base),offset=k%per_base;
                    temp=base_counter_size;
                }
                else{
                    word_index=(k/per_magnify),offset=k%per_magnify;
                    temp=magnify_counter_size;
                }

                if((allcounter[i][j][word_index] >> (offset*temp)) & (uint64)(1<<(temp-1))){
                    used++;
                }
            }
            rate[j]+=(double)used/(dep*wid);
            cou[j]++;
        }
    }

    for(int i=0;i<15;i++){
        if(cou[i]){
            cout<<"the overflow rate of level "<<i<<": "<<rate[i]/cou[i]<<endl;
        }
    }
}

//calculate the use rate of each level
void PCMSketch::useRate(){
    int cou[15];
    double rate[15];
    for(int i=0;i<15;i++){
        cou[i]=0;
        rate[i]=0.0;
    }
    int dep,wid,used;
    for(int i=0;i<d;i++){
        for(int j=0;j<allcounter[i].size();j++){
            dep=depth[i][j];wid=width[i][j];
            used=0;
            for(int k=0;k<dep*wid;k++){
                int word_index,offset,temp;
                if(j==0){
                    word_index=(k/per_base),offset=k%per_base;
                    temp=base_counter_size;
                }
                else{
                    word_index=(k/per_magnify),offset=k%per_magnify;
                    temp=magnify_counter_size;
                }
                if((allcounter[i][j][word_index] >> (temp*offset)) & (uint64)((1<<temp)-1)){
                    used++;
                }
            }
            rate[j]+=(double)used/(dep*wid);
            cou[j]++;
        }
    }
    for(int i=0;i<15;i++){
        if(cou[i]){
            cout<<"the userate of level "<<i<<": "<<rate[i]/cou[i]<<endl;
        }
    }
}
//get the physical location
int PCMSketch::getrealpos(int hash_value1,int hash_value2,int level,int num){
    //calculate the area
    int area,last=level-1,realpos;
    bool up=((hash_value1<<1)<depth[num][last]),left=((hash_value2<<1)<width[num][last]);
    if(up&&left){
        area=0;
    }
    else if(up&&!left){
        area=1;
    }
    else if(!up&&left){
        area=2;
    }
    else{
        area=3;
    }
    switch(area){
        case 0 :
        {
            realpos=hash_value1*width[num][level]+hash_value2;
            break;
        }
        case 1 :
        {
            realpos=depth[num][level]*width[num][level]+hash_value1*(width[num][last]-width[num][level])+(hash_value2-width[num][level]);
            break;
        }
        case 2 :
        {
            realpos=depth[num][level]*width[num][level]+depth[num][level]*(width[num][last]-width[num][level])+(hash_value1-depth[num][level])*width[num][level]+hash_value2;
            break;
        }
        default :
        {
            realpos=depth[num][level]*width[num][level]+depth[num][level]*(width[num][last]-width[num][level])+(depth[num][last]-depth[num][level])*width[num][level]+
            (hash_value1-depth[num][level])*(width[num][last]-width[num][level])+(hash_value2-width[num][level]);
        }
    }
    return realpos;
}
void PCMSketch::get_down_pos(int pos, int level,int num,vector<int>& res){//from level to level-1
    int x=pos/width[num][level],y=pos%width[num][level];

    int x1=x,y1=y;
    if(x1<depth[num][level-1]&&y1<width[num][level-1])res.push_back(x1*width[num][level-1]+y1);

    x1=(depth[num][level]-1-x),y1=((width[num][level]<<1)-1-y);
    if(x1<depth[num][level-1]&&y1<width[num][level-1])res.push_back(x1*width[num][level-1]+y1);

    x1=((depth[num][level]<<1)-1-x);y1=(width[num][level]-1-y);
    if(x1<depth[num][level-1]&&y1<width[num][level-1])res.push_back(x1*width[num][level-1]+y1);

    x1=depth[num][level]+x;y1=width[num][level]+y;
    if(x1<depth[num][level-1]&&y1<width[num][level-1])res.push_back(x1*width[num][level-1]+y1);
}
int PCMSketch::getpos(int pos,int level,int num,int& area){//from level-1 to level
//    string str=to_string(pos);
//    return (bobhash[num]->run(str.c_str(), strlen(str.c_str())))%(depth[num][level]*width[num][level]);
//
//    int x=pos/width[num][level-1],y=pos%width[num][level-1];
//    return (x>>1)*width[num][level]+(y>>1);

    int x=pos/width[num][level-1],y=pos%width[num][level-1];
    //calculate the area
    bool up=x<depth[num][level],left=y<width[num][level];
    int newx,newy;
    if(up&&left){
        area=0;
        newx=x;newy=y;
    }
    else if(up&&!left){
        area=1;
        newx=depth[num][level]-1-x;
        newy=(width[num][level]<<1)-1-y;
    }
    else if(!up&&left){
        area=2;
        newx=(depth[num][level]<<1)-1-x;
        newy=width[num][level]-1-y;
    }
    else{
        area=3;
        newx=x-depth[num][level];
        newy=y-width[num][level];
    }
    return newx*width[num][level]+newy;
//      return (pos>>2);
}
//w means the number of counter and _d means the number of hash function
PCMSketch::PCMSketch(int w, int _d,int _base_counter_size,int _magnify_counter_size,double up,double buttom)
{
	d = _d;
	this->w=w;
	wordNum=0;
	leftRatio=0;
	appendDuration=0;
	appendTime=0;
	up_bound=(int)(up*w);
	buttom_bound=(int)(buttom*w);
	//cout<<buttom_bound<<" "<<up_bound<<endl;

	base_counter_size=_base_counter_size;
	magnify_counter_size=_magnify_counter_size;
	per_base=64/base_counter_size;
	per_magnify=64/magnify_counter_size;
	max_base=((1<<(base_counter_size-1))-1);
	max_magnify=((1<<(magnify_counter_size-1))-1);

	overflow_count=new int[d];
	base_level=new int[d];
	memset(base_level,0,sizeof(int)*d);
	memset(overflow_count,0,sizeof(int)*d);
	//get the shape of the hash matrix
	for(int i=0;i<d;i++){
        vector<int> mid1,mid2;
        depth.push_back(mid1);
        width.push_back(mid2);
	}
	depth[0].push_back(sqrt(w));
	width[0].push_back(w/depth[0][0]);
    for(int i=1;i<=(d>>1);i++){
        int divd=1<<(2*i);
        int temp1=sqrt(w)/divd;
        if(temp1==0)temp1=1;
        int temp2=w/temp1;
        width[2*i-1].push_back(temp1);
        depth[2*i-1].push_back(temp2);

        width[2*i].push_back(temp2);
        depth[2*i].push_back(temp1);
    }

    for(int j=0;j<d;j++){
        int dep=depth[j][0],wid=width[j][0],counter_num;
        vector<uint64*> counter;
        vector<bool*> eachFlag;
        bool next=true;
        for(int i = 0; i < 15; i++){
            if(!next)break;
            counter_num= dep*wid;
            int word_num;
	        if(i==0){
                word_num=(counter_num/per_base)+(int)(counter_num%per_base!=0);
                wordNum=max((counter_num/64)+(int)(counter_num%64!=0),wordNum);
	        }
	        else word_num=(counter_num/per_magnify)+(int)(counter_num%per_magnify!=0);
            if(i!=0){
                depth[j].push_back(dep);
                width[j].push_back(wid);
            }

		    uint64* mid = new uint64[word_num];
		    memset(mid, 0, sizeof(uint64) * word_num);
		    if(counter_num==1)next=false;


		    dep=((dep+1)>>1);
            wid=((wid+1)>>1);
		    counter.push_back(mid);
		    if(i<PARALLELCOUNTER-1){
                bool *flagMid=new bool[counter_num];
                memset(flagMid,false,sizeof(bool)*counter_num);
                eachFlag.push_back(flagMid);
		    }
	   }
	   allcounter.push_back(counter);
	   accessedFlag.push_back(eachFlag);
    }
	for(int i = 0; i < 15; i++)
	{
		bobhash[i] = new BOBHash64(i + 1000);
	}
}

void PCMSketch::Insert(string str1,string str2,int weight,unsigned long& hashTime,unsigned long& updateTime)
{
    uint64 hash_value1,hash_value2,value;

    uint64 hash1,hash2,helper=(((uint64)0x1)<<32)-1,useless;
    unsigned long long h1,h2;
    unsigned int hashval1[16],hashval2[16];

    int i=0;
    timeval t_start, t_end;
    gettimeofday( &t_start, NULL);
    while(true){
        //for(int i=0;i<3;i++){
            h1=(bobhash[i]->run(str1.c_str(), strlen(str1.c_str())));
            h2=(bobhash[i]->run(str2.c_str(), strlen(str2.c_str())));
        //}
        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        h1>>=32;h2>>=32;i++;

        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        i++;
        if(i>=d)break;
    }
    gettimeofday( &t_end, NULL);
    hashTime+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;

    int pos_all[d];
    for(int i=0;i<d;i++){
        hash_value1 = hashval1[i]%depth[i][0];
	    hash_value2 = hashval2[i]%width[i][0];
	    pos_all[i]=hash_value1*width[i][0]+hash_value2;
    }
    gettimeofday( &t_start, NULL);
    for(int i=0;i<d;i++){
        int pos=pos_all[i];
	    int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

	    uint64 base=0;
	    int word_index=(pos/64),offset=(pos%64);
	    for(int j=0;j<base_level[i];j++){
            base+=(((allcounter[i][j][word_index] >> offset) & (0x1))<<j);
	    }
	    base+=weight;
	    int temp=(1<<base_level[i]);
	    int overflow=base/temp,left=base%temp;
	    for(int j=0;j<base_level[i];j++){
            allcounter[i][j][word_index]&=(~((uint64)0x1 << offset));
            allcounter[i][j][word_index]|= ((uint64)(left&0x1) << offset);
            left>>=1;
	    }

        value = (allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & max_base;

        int now=overflow+value;
        overflow=(now/(max_base+1)),left=now%(max_base+1);

        bool overflow_flag=((allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & (uint64)(1<<(base_counter_size-1)));
        allcounter[i][base_level[i]][my_word_index]&=(~((uint64)max_base << (counter_offset*base_counter_size)));
        allcounter[i][base_level[i]][my_word_index]|= ((uint64)left << (counter_offset*base_counter_size));
        if (overflow!=0)
        {
            allcounter[i][base_level[i]][my_word_index]|= ((uint64)(max_base+1) << (counter_offset*base_counter_size));
            carry(overflow,i,pos);
        }
        if(!overflow_flag&&((allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & (uint64)(max_base+1)))
                overflow_count[i]++;
        if(overflow_count[i]>up_bound){
            down_flow(i);
        }
    }
    gettimeofday( &t_end, NULL);
    updateTime+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;
}

void PCMSketch::Delete(string str1,string str2,int weight){
    uint64 hash_value1,hash_value2,value;

    uint64 hash1,hash2,helper=(((uint64)0x1)<<32)-1,useless;
    unsigned long long h1,h2;
    unsigned int hashval1[16],hashval2[16];
    int i=0;
    while(true){
        h1=(bobhash[i]->run(str1.c_str(), strlen(str1.c_str())));
        h2=(bobhash[i]->run(str2.c_str(), strlen(str2.c_str())));
        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        h1>>=32;h2>>=32;i++;

        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        i++;
        if(i>=d)break;
    }
    int min_value=1<<30;
    int pos_all[d];
    for(int i=0;i<d;i++){
        hash_value1 = hashval1[i]%depth[i][0];
	    hash_value2 = hashval2[i]%width[i][0];
	    pos_all[i]=hash_value1*width[i][0]+hash_value2;
    }
    for(int i=0;i<d;i++){
        int pos=pos_all[i];

	    int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

        value = (allcounter[i][0][my_word_index] >> (counter_offset*base_counter_size)) & max_base;
        int overflow=(weight/(max_base+1)),left=(weight%(max_base+1));
        bool overflow_flag=((allcounter[i][0][my_word_index] >> (counter_offset*base_counter_size)) & (uint64)(1<<(base_counter_size-1)));

        if(left>value){
            //not enough to sub and no overflow
            if(!overflow_flag){
                allcounter[i][0][my_word_index]&=(~((uint64)max_base << (counter_offset*base_counter_size)));
                //cout<<"error\n";
                return;
            }
            overflow++;
            value+=(max_base+1);
        }
        int now=value-left;
        allcounter[i][0][my_word_index]&=(~((uint64)max_base << (counter_offset*base_counter_size)));
        allcounter[i][0][my_word_index]|= ((uint64)now << (counter_offset*base_counter_size));

        if (overflow!=0)
        {
            if(!overflow_flag){
                allcounter[i][0][my_word_index]&=(~((uint64)max_base << (counter_offset*base_counter_size)));
                //cout<<"error\n";
                return;
            }
            carryDelete(overflow,i,pos,1);

        }
    }
}
int PCMSketch::nodeQuery(string str1, bool type){
    unsigned int m=1<<30;

    uint64 hash_value1,hash_value2,value;

    uint64 hash1,hash2,helper=(((uint64)0x1)<<32)-1,useless;
    unsigned long long h1,h2;
    unsigned int hashval1[16],hashval2[16];



    int i=0;
    while(true){
        h1=(bobhash[i]->run(str1.c_str(), strlen(str1.c_str())));
        hashval1[i]=h1&(helper);
        h1>>=32;i++;

        hashval1[i]=h1&(helper);
        i++;
        if(i>=d)break;
    }
	for (int i = 0; i < d; i++){
        unsigned int inner=0;
		if(type){
            if(i!=0&&i%2==0)continue;
            hash_value1 = hashval1[i]%depth[i][0];
            for(int j=0;j<width[i][0];j++){
                int pos=hash_value1*width[i][0]+j;
//                pos=getrealpos(hash_value1,j,1,i);
	            int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

                int word_index=(pos/64),offset=(pos%64);
	            for(int j=0;j<base_level[i];j++){
                    inner+=(((allcounter[i][j][word_index] >> offset) & (0x1))<<j);
	            }

                inner+=(((allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & max_base)<<base_level[i]);

                if ((allcounter[i][base_level[i]][my_word_index]>>((counter_offset*base_counter_size)+base_counter_size-1))&0x1){
                    inner+=get_value(i,pos);
                }
            }
		}
		else{
            if(i%2)continue;
		    hash_value1 = hashval1[i]%width[i][0];
            for(int j=0;j<depth[i][0];j++){

	           int pos=j*width[i][0]+hash_value1;
//	           pos=getrealpos(j,hash_value1,1,i);
	           int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

               int word_index=(pos/64),offset=(pos%64);
               for(int j=0;j<base_level[i];j++){
                    inner+=(((allcounter[i][j][word_index] >> offset) & (0x1))<<j);
               }

               inner+=(((allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & max_base)<<base_level[i]);
               if ((allcounter[i][base_level[i]][my_word_index]>>((counter_offset*base_counter_size)+base_counter_size-1))&0x1){
                   inner+=get_value(i,pos);
               }
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

int PCMSketch::edgeQuery(string str1,string str2)
{
	int min_value = 1 << 30;

    uint64 hash_value1,hash_value2,value;

    uint64 hash1,hash2,helper=(((uint64)0x1)<<32)-1,useless;
    unsigned long long h1,h2;
    unsigned int hashval1[16],hashval2[16];

    int i=0;
    while(true){
        h1=(bobhash[i]->run(str1.c_str(), strlen(str1.c_str())));
        h2=(bobhash[i]->run(str2.c_str(), strlen(str2.c_str())));
        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        h1>>=32;h2>>=32;i++;

        hashval1[i]=h1&(helper);hashval2[i]=h2&helper;
        i++;
        if(i>=d)break;
    }

    for(int i=0;i<d;i++){
        hash_value1 = hashval1[i]%depth[i][0];
	    hash_value2 = hashval2[i]%width[i][0];
	    int pos=hash_value1*width[i][0]+hash_value2;
//	    pos=getrealpos(hash_value1,hash_value2,1+base_level[i],i);
	    int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

	    int word_index=(pos/64),offset=(pos%64);
	    value=0;
	    for(int j=0;j<base_level[i];j++){
            value+=(((allcounter[i][j][word_index] >> offset) & (0x1))<<j);
	    }

        value+=(((allcounter[i][base_level[i]][my_word_index] >> (counter_offset*base_counter_size)) & max_base)<<base_level[i]);

        if ((allcounter[i][base_level[i]][my_word_index]>>((counter_offset*base_counter_size)+base_counter_size-1))&0x1)
        {
            value+=get_value(i,pos);
        }
        if(value<min_value)min_value=value;
        if(min_value==0)return 0;
    }
	return min_value;

}

void PCMSketch::resetOverflowFlag(int pos,int level,int num){
    vector<int> posDown;

    if(level>base_level[num])get_down_pos(pos,level,num,posDown);
    else posDown.push_back(pos);

    int per,counter_size,my_word_index,counter_offset,max_here;
    if(level==1){per=per_base;counter_size=base_counter_size;max_here=max_base;}
    else{per=per_magnify;counter_size=magnify_counter_size;max_here=max_magnify;}

    for(int z=0;z<posDown.size();z++){
        int mid=posDown[z];
        my_word_index=mid/per,counter_offset=mid%per;
        if(((allcounter[num][level-1][my_word_index] >> (counter_offset*counter_size)) & (uint64)(1<<(counter_size-1)))){
            if(level-1==base_level[num])overflow_count[num]--;
            allcounter[num][level-1][my_word_index]&= (~((uint64)(max_here+1) << (counter_offset*counter_size)));
            if(!((allcounter[num][level-1][my_word_index] >> (counter_size*counter_offset)) & (uint64)((1<<counter_size)-1))&&level-1>0)
                resetOverflowFlag(mid,level-1,num);
        }
    }
}
void PCMSketch::carryDelete(int overflow,int num, int pos,int level){
    int area;
    if(level>base_level[num]){
        pos=getpos(pos,level,num,area);
    }
    int word_index=(pos/per_magnify),offset=(pos%per_magnify);

    int value = (allcounter[num][level][word_index] >> (offset*magnify_counter_size)) & max_magnify;
    int temp=overflow;
    overflow=(temp/(max_magnify+1));
    int left=temp%(max_magnify+1);
    bool overflow_flag=((allcounter[num][level][word_index] >> (offset*magnify_counter_size)) & (uint64)(1<<(magnify_counter_size-1)));
    bool processed=false;
    if(left>value){
        //not enough to sub and no overflow
        if(!overflow_flag){
            allcounter[num][level][word_index]&=(~((uint64)max_magnify << (offset*magnify_counter_size)));
            //cout<<"error\n";
            processed=true;
        }
        else{
            overflow++;
            value+=(max_magnify+1);
        }
    }
    if(!processed){
        int now=value-left;

        allcounter[num][level][word_index]&=(~((uint64)max_magnify << (offset*magnify_counter_size)));
        allcounter[num][level][word_index]|= ((uint64)now << (offset*magnify_counter_size));

        if(overflow>0){
            if(!overflow_flag){//enough to sub this level but not enough to sub more
                allcounter[num][level][word_index]&=(~((uint64)max_magnify << (offset*magnify_counter_size)));
                //cout<<"error\n";
            }
            else carryDelete(overflow,num,pos,level+1);
        }
    }

    if(!((allcounter[num][level][word_index] >> (offset*magnify_counter_size))& (uint64)((1<<magnify_counter_size)-1))){
        resetOverflowFlag(pos,level,num);
//        if(overflow_count[num]<buttom_bound&&base_level[num]>0)up_flow(num);
    }
}


void PCMSketch::carry(int overflow,int num,int pos)
{
    int area;
	for(int i = base_level[num]+1; i < allcounter[num].size(); i++)
	{
	    pos=getpos(pos,i,num,area);
	    int word_index=(pos/per_magnify),offset=(pos%per_magnify);

		int value = (allcounter[num][i][word_index] >> (offset*magnify_counter_size)) & max_magnify;
		int now=value+overflow,overflow=(now/(max_magnify+1)),left=now%(max_magnify+1);
		bool overflow_flag=((allcounter[num][i][word_index] >> (offset*magnify_counter_size)) & (uint64)(1<<(magnify_counter_size-1)));

		allcounter[num][i][word_index]&=(~((uint64)max_magnify << (offset*magnify_counter_size)));
		allcounter[num][i][word_index]|= ((uint64)left << (offset*magnify_counter_size));

		if(overflow==0){
            return;
		}
		allcounter[num][i][word_index]|= ((uint64)(max_magnify+1) << (offset*magnify_counter_size));
	}
}

int PCMSketch::get_value(int num,int pos)
{
    int res=0,area;
    for(int i = base_level[num]+1; i < allcounter[num].size(); i++)
	{
	    pos=getpos(pos,i,num,area);
	    int word_index=(pos/per_magnify),offset=(pos%per_magnify);
		int value = (allcounter[num][i][word_index] >> (offset*magnify_counter_size)) & max_magnify;

		res+=value<<(base_counter_size-1+(i-base_level[num]-1)*(magnify_counter_size-1)+base_level[num]);
		if((allcounter[num][i][word_index]>>((offset*magnify_counter_size)+magnify_counter_size-1))&0x1){
            continue;
		}
		break;
	}
	return res;
}

#endif //_PCMSKETCH_H
