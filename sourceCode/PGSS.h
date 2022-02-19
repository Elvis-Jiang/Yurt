#ifndef PGSS_H_INCLUDED
#define PGSS_H_INCLUDED
#include "querysupportstruct.h"
using namespace std;

PGSS::PGSS(int _width, int range, int p_num, int size,int f_num, int _base_counter_size, int _magnify_counter_size, double _up_bound, int _flag_size, bool usehashtable, int TableSize)//the side length of matrix, the length of hash addtress list, the number of candidate bucekt
// the number of rooms, whether to use hash table, and the size of the table.
// Hash table which stores the original nodes can be omitted if not needed. For nodequery,
// reachability, edgequery not needed. But needed for triangel counting, degree query, and successor / precursor queries.
{
    bufferSize=0;
    appendDuration=0;
    appendTime=0;

    base_counter_size=_base_counter_size;
    flag_size=_flag_size;
    max_flag=((1<<flag_size)-1);
	magnify_counter_size=_magnify_counter_size;
	per_base=64/base_counter_size;
	per_magnify=64/magnify_counter_size;
	max_base=((1<<(base_counter_size-1))-1);
	max_magnify=((1<<(magnify_counter_size-flag_size))-1);
	max_value_base=((1<<base_counter_size)-1);
	max_value_magnify=((1<<magnify_counter_size)-1);
    base_level=0;overflow_count=0;

	w = _width;
	r = range; /* r x r mapped baskets */
	p = p_num; /*candidate buckets*/
	s = size; /*multiple rooms*/
	f = f_num; /*finger print lenth*/
	n = 0;
	up_bound=(int)(_up_bound*w*w);
	value = new Pbasket[w*w];
	useT=usehashtable;
	tablesize=TableSize;
	memset(value, 0, sizeof(Pbasket)*w*w);
	collision=new int[15];
	memset(collision, 0, sizeof(int)*15);
	wordNum=((w*w)/64)+(int)((w*w)%64!=0);


	int wid=w,counter_num;
    bool next=true;
	for(int i = 0; i < 15; i++){
        if(!next)break;
        counter_num= wid*wid;
        int word_num;
        if(i==0){
            word_num=(counter_num/per_base)+(int)(counter_num%per_base!=0);
        }
        else word_num=(counter_num/per_magnify)+(int)(counter_num%per_magnify!=0);
        width.push_back(wid);

        uint64* mid = new uint64[word_num];
        memset(mid, 0, sizeof(uint64) * word_num);
        if(counter_num==1)next=false;

        wid=((wid+1)>>1);
        allcounter.push_back(mid);
        if(i<PARALLELCOUNTER-1){
            bool *flagMid=new bool[counter_num];
            memset(flagMid,false,sizeof(bool)*counter_num);
            accessedFlag.push_back(flagMid);
        }
    }

	if(usehashtable)
	mapTable.init(tablesize);
}
void PGSS::show_info(){
    cout<<"the overflow ratio of the first layer is "<<(double)(overflow_count)/(w*w)<<endl;
    for(int i=1;i<allcounter.size();i++){
        int num=width[i]*width[i];
        double count=0;
        for(int j=0;j<num;j++){
           int word_index=(j/per_magnify),offset=(j%per_magnify);
           if((allcounter[i][word_index] >> (offset*magnify_counter_size))&12)count+=1;
        }
        cout<<"the overflow ration of the "<<i<<" layer is "<<count/num<<endl;
    }
    cout<<endl;
    for(int i=0;i<allcounter.size();i++){
        cout<<"the collision ration of the "<<i<<" layer is "<<collision[i]<<endl;
    }
}
void PGSS::matrixUserate(){
    double used=0;
    for(int i=0;i<w*w;i++){
        for(int j=0;j<s;j++){
            if(value[i].src[j])used+=1;
        }
    }
    cout<<bufferSize<<" "<<used/(w*w*s)<<endl;
}
bool PGSS::occupied(int pos){
    pos=getpos(pos,base_level+1);
    int word_index=(pos/per_magnify),offset=(pos%per_magnify);
    return ((allcounter[base_level+1][word_index] >> (offset*magnify_counter_size))&max_value_magnify);
}
//append the structure
void PGSS::down_flow(){
    cout<<"down flow\n";
    timeval t_start, t_end;
    gettimeofday( &t_start, NULL);
    appendTime+=1;
    uint64 *a=new uint64[wordNum];

    memset(a,0,sizeof(uint64)*wordNum);

    allcounter.insert(allcounter.begin()+base_level,a);
    int wid=width[0];
    width.insert(width.begin()+base_level,wid);

    base_level++;

    for(int i=base_level-1;i<allcounter.size();i++){//down flow
        //the base level, we pull directly
        int counter_num=width[i]*width[i],level=i;
        if(level==base_level-1){
            for(int pos=0;pos<counter_num;pos++){
                int my_word_index=(pos/per_base),counter_offset=(pos%per_base);
                uint64 downValue=((allcounter[level+1][my_word_index] >> (counter_offset*base_counter_size)) & max_base);

                int word_index=(pos/64),offset=(pos%64);
                uint64 left=(downValue>>1),value=(downValue&0x1);

                allcounter[level][word_index]|= (value << offset);
            }
        }
        //the mid level, we move, check overflow and pull
        else{
            int per,counter_size,max_here;
            if(level==base_level){
                per=per_base;
                counter_size=base_counter_size;
                max_here=max_base;
            }
            else{
                per=per_magnify;
                counter_size=magnify_counter_size;
                max_here=max_magnify;
            }

            for(int pos=0;pos<counter_num;pos++){
                int my_word_index=(pos/per),counter_offset=(pos%per);
                //we get the value and right move
                uint64 valueNow=(allcounter[level][my_word_index] >> (counter_offset*counter_size)) & max_here;
                valueNow>>=1;
                //we get the overflow flag
                int overflow_flag;
                if(level==base_level){
                    overflow_flag=((allcounter[level][my_word_index] >> (counter_offset*counter_size+counter_size-1)) &0x1);
                }
                else{
                    overflow_flag=((allcounter[level][my_word_index] >> (counter_offset*counter_size+counter_size-flag_size)) &max_flag);
                }
                //we check the overflow and pull the res from the upper level
                if(overflow_flag>0){
                    int posUpper=(getpos(pos,level+1)+overflow_flag-1)%(width[level+1]*width[level+1]);
                    //we need to check the upper's situation and make sure we have to overflow
                    int upper_word_index=(posUpper/per_magnify),upper_counter_offset=(posUpper%per_magnify);
                    int upper_overflow_flag=((allcounter[level+1][upper_word_index] >>
                                          (upper_counter_offset*magnify_counter_size+magnify_counter_size-flag_size)) &max_flag);
                    uint64 upper_value=(allcounter[level+1][upper_word_index] >>
                                    (upper_counter_offset*magnify_counter_size)) & max_magnify;

                    uint64 left=(upper_value>>1),down_value=(upper_value&0x1);
                    bool upper_exist=(left||upper_overflow_flag>0);
                    //down flow the bit and maybe we got to reset the flag
                    if(level==base_level){
                        valueNow|=(down_value<<(counter_size-2));
                        if(!upper_exist){
                            allcounter[level][my_word_index]&= (~((uint64)(max_here+1) << (counter_offset*counter_size)));
                            overflow_count--;
                        }
                    }
                    else{
                        valueNow|=(down_value<<(counter_size-flag_size-1));
                        if(!upper_exist){
                            allcounter[level][my_word_index]&= (~((uint64)(max_flag) << (counter_offset*counter_size+counter_size-flag_size)));
                        }
                    }
                }

                //we reset the value
                allcounter[level][my_word_index]&=(~((uint64)max_here << (counter_offset*counter_size)));
                allcounter[level][my_word_index]|= ((uint64)valueNow << (counter_offset*counter_size));
            }
        }
    }

    gettimeofday( &t_end, NULL);
    appendDuration+=(t_end.tv_sec-t_start.tv_sec)*1000.0 +
                    (t_end.tv_usec-t_start.tv_usec)/1000.0;
}
void PGSS::get_down_pos(int pos, int level,vector<int>& res){//from level to level-1
    int x=pos/width[level],y=pos%width[level];

    int x1=x,y1=y;
    if(x1<width[level-1]&&y1<width[level-1])res.push_back(x1*width[level-1]+y1);

    x1=(width[level]-1-x),y1=((width[level]<<1)-1-y);
    if(x1<width[level-1]&&y1<width[level-1])res.push_back(x1*width[level-1]+y1);

    x1=((width[level]<<1)-1-x);y1=(width[level]-1-y);
    if(x1<width[level-1]&&y1<width[level-1])res.push_back(x1*width[level-1]+y1);

    x1=width[level]+x;y1=width[level]+y;
    if(x1<width[level-1]&&y1<width[level-1])res.push_back(x1*width[level-1]+y1);
}
int PGSS::getpos(int pos,int level){//from level-1 to level
//    string str=to_string(pos);
//    return (bobhash[num]->run(str.c_str(), strlen(str.c_str())))%(depth[num][level]*width[num][level]);
//
//    int x=pos/width[num][level-1],y=pos%width[num][level-1];
//    return (x>>1)*width[num][level]+(y>>1);

    int x=pos/width[level-1],y=pos%width[level-1];
    //calculate the area
    bool up=x<width[level],left=y<width[level];
    int newx,newy;
    if(up&&left){
        newx=x;newy=y;
    }
    else if(up&&!left){
        newx=width[level]-1-x;
        newy=(width[level]<<1)-1-y;
    }
    else if(!up&&left){
        newx=(width[level]<<1)-1-x;
        newy=width[level]-1-y;
    }
    else{
        newx=x-width[level];
        newy=y-width[level];
    }
    return newx*width[level]+newy;
//      return (pos>>2);
}
//insert into the second and high layer of the pyramid
void PGSS::carry(int overflow,int pos)
{
    //insert into the upper layers
    int last=0;
    int word_index,offset;
	for(int i = base_level+1; i < allcounter.size(); i++)
	{
	    //set the overflow flag of the last layer
	    if(last){
            allcounter[i-1][word_index]|= ((uint64)(last) << (offset*magnify_counter_size+magnify_counter_size-flag_size));
	    }
	    word_index=(pos/per_magnify);offset=(pos%per_magnify);

		int value = (allcounter[i][word_index] >> (offset*magnify_counter_size)) & max_magnify;
		int now=value+overflow,overflow=(now/(max_magnify+1)),left=now%(max_magnify+1);
		//if overflowed already, we get the overflow position
		int overflow_flag=((allcounter[i][word_index] >> (offset*magnify_counter_size+magnify_counter_size-flag_size)) &max_flag);

		allcounter[i][word_index]&=(~((uint64)max_magnify << (offset*magnify_counter_size)));
		allcounter[i][word_index]|= ((uint64)left << (offset*magnify_counter_size));

		if(overflow==0){
            return;
		}
		//select the best overflow position
		pos=getpos(pos,i+1);
		int lastpos=pos;
		last=-1;
		int max_size=width[i+1]*width[i+1];
		if(overflow_flag){
            last=overflow_flag;
            pos=(pos+last-1)%(max_size);
            continue;
		}
		//if it is the first time to overflow, we select a good position
		for(int j=0;j<max_flag;j++){
		    int posnow=(lastpos+j)%max_size;
		    int idx=(posnow/per_magnify),off=(posnow%per_magnify);
		    //cout<<((allcounter[i+1][idx] >> (off*magnify_counter_size))&max_value_magnify)<<endl;
            if(!((allcounter[i+1][idx] >> (off*magnify_counter_size))&max_value_magnify)){
                pos=posnow;
                last=j+1;
                break;
            }
		}
		if(last==-1){
            collision[i-base_level]+=1;
            last=(rand()%max_flag)+1;
            lastpos+=(last-1);
		}
	}
}
//insert into the first layer of the pyramid
void PGSS::insertMatrix(int pos,int weight){
    uint64 base=0;
    int my_word_index=(pos/per_base),counter_offset=(pos%per_base);
    int word_index=(pos/64),offset=(pos%64);
    for(int j=0;j<base_level;j++){
        base+=(((allcounter[j][word_index] >> offset) & (0x1))<<j);
    }
    base+=weight;
    int temp=(1<<base_level);
    int overflow=base/temp,left=base%temp;
    for(int j=0;j<base_level;j++){
        allcounter[j][word_index]&=(~((uint64)0x1 << offset));
        allcounter[j][word_index]|= ((uint64)(left&0x1) << offset);
        left>>=1;
    }

    int val = (allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & max_base;


    //now, we insert into the pyramid
    int now=overflow+val;
    overflow=(now/(max_base+1)),left=now%(max_base+1);

    bool overflow_flag=((allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & (uint64)(1<<(base_counter_size-1)));
    allcounter[base_level][my_word_index]&=(~((uint64)max_base << (counter_offset*base_counter_size)));
    allcounter[base_level][my_word_index]|= ((uint64)left << (counter_offset*base_counter_size));

    if (overflow!=0)
    {
        allcounter[base_level][my_word_index]|= ((uint64)(max_base+1) << (counter_offset*base_counter_size));
        //the match strategy of the first level to the second is fixed
        pos=getpos(pos,base_level+1);
        carry(overflow,pos);
    }

    if(!overflow_flag&&((allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & (uint64)(max_base+1)))
            overflow_count++;
    if(overflow_count>up_bound){
        down_flow();
    }
}
void PGSS::cleanupBuffer()
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
void PGSS::insert(string s1, string s2,int weight,uint64& insertTime)// s1 is the ID of the source node, s2 is the ID of the destination node, weight is the edge weight.
{
        timeval t_start, t_end;
        gettimeofday( &t_start, NULL);
		unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
		unsigned int hash2 = (*hfunc[0])((unsigned char*)(s2.c_str()), s2.length());
		unsigned int tmp = pow(2,f)-1;

		unsigned short g1 = hash1 & tmp;
		if(g1==0) g1+=1;
		unsigned int h1 = (hash1>>f)%w;
		unsigned short g2 = hash2 & tmp;
		if(g2==0) g2+=1;
		unsigned int h2 = (hash2>>f)%w;

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
		bool inserted=false,findpos=false;
		long key = g1+g2;
		int emptyIdx=-1,emptyNotOverflowIdx=-1,itempos=-1;
		for(int i=0;i<p;i++)
		{
			key = (key*timer+prime)%bigger_p;
			int index = key%(r*r);
			int index1 = index/r;
			int index2 = index%r;
			int p1 = (h1+tmp1[index1])%w;
			int p2 = (h2+tmp2[index2])%w;


			int pos = p1*w + p2;
			int my_word_index=(pos/per_base),counter_offset=(pos%per_base);
			for (int j = 0; j < s; j++)
			{
			    //find a empty bucket and we record
			    if (value[pos].src[j] == 0)
				{
				    if(emptyIdx==-1)
                        emptyIdx=index;
					if(emptyNotOverflowIdx==-1&&!occupied(pos))
					    emptyNotOverflowIdx=index;
				}
				else if ( value[pos].idx[j] == (index2 | (index1 << 8)) && (value[pos].src[j]== g1) && (value[pos].dst[j] == g2) )
				{
				    itempos=pos;
				    //we do not need to move if we have already overflow or we do not need to overflow
					if(((allcounter[base_level][my_word_index]>>((counter_offset*base_counter_size)+base_counter_size-1))&0x1)||
                         (((allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & max_base)+weight)<max_base){
                         insertMatrix(pos,weight);
					     inserted = true;
                         break;
                    }
                    //we got to move and we have find a good position
                    else{
                        if(emptyNotOverflowIdx!=-1){
                            findpos=true;
                            break;
                        }
                    }
				}
			}
			if(inserted||findpos)
				break;
		}
		//we find the item but it need a new position
		if(!inserted&&itempos!=-1){
            int my_word_index=(itempos/per_base),counter_offset=(itempos%per_base);



            //we find a not overflow position
            if(emptyNotOverflowIdx!=-1){
                inserted = true;

                //clear the old position
                value[itempos].src[0]=0;
                value[itempos].dst[0]=0;
                value[itempos].idx[0]=0;
                weight+=((allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & max_base);
                allcounter[base_level][my_word_index]&=(~((uint64)max_value_base << (counter_offset*base_counter_size)));

                //get the new position
                int index1 = emptyNotOverflowIdx/r;
			    int index2 = emptyNotOverflowIdx%r;
			    int p1 = (h1+tmp1[index1])%w;
			    int p2 = (h2+tmp2[index2])%w;
			    int pos = p1*w + p2;
			    //insert into the new position
			    value[pos].src[0]=g1;
			    value[pos].dst[0]=g2;
			    value[pos].idx[0]=((index1<<8)|index2);
			    insertMatrix(pos,weight);
            }
            //we do not find a good one, so we have to keep it
            else{
                insertMatrix(itempos,weight);
                //count the collision number of the counter
                collision[0]++;
                inserted=true;
            }
		}
		linknode *frontNode=NULL;
		//we do not find the item in the matrix and not empty bucket, which means the edge maybe in the buffer
		if(!inserted)
		{
			map<unsigned int, int>::iterator it = index.find(k1);
			if(it!=index.end())
			{
				int tag = it->second;
				linknode* node = buffer[tag];
				while(node)
				{
				    frontNode=node;
				    //we find the item in the buffer
					if (node->key == k2)
					{
						node->weight += weight;
						inserted=true;
						break;
					}
					node = node->next;
				}
			}
		}
		//this is a new item and we try to insert it into the matrix
		if(!inserted){
		    int insertIdx;
		    if(weight>max_base){
                if(emptyNotOverflowIdx!=-1){
                    inserted=true;
                    insertIdx=emptyNotOverflowIdx;
                }
                else if(emptyIdx!=-1){
                    inserted=true;
                    insertIdx=emptyIdx;
                }
		    }
		    else{
                if(emptyIdx!=-1){
                    inserted=true;
                    insertIdx=emptyIdx;
                }
		    }
            if(inserted){
                int index1 = insertIdx/r;
                int index2 = insertIdx%r;
                int p1 = (h1+tmp1[index1])%w;
                int p2 = (h2+tmp2[index2])%w;
                int pos = p1*w + p2;
                //insert into the new position
                value[pos].src[0]=g1;
                value[pos].dst[0]=g2;
                value[pos].idx[0]=((index1<<8)|index2);
                insertMatrix(pos,weight);
            }
		}
		//we have to insert the item into the buffer
		if(!inserted){
            //we just insert into the existed linked list
            if(frontNode){
                linknode* ins = new linknode;
                bufferSize++;
                ins->key = k2;
                ins->weight = weight;
                ins->next = NULL;
                frontNode->next = ins;
            }
            //we got to make a new linked list
            else{
                index[k1] = n;
				n++;
				linknode* node = new linknode;
				bufferSize++;
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
int PGSS::get_value(int pos)
{
    int res=0;
    for(int i = base_level+1; i < allcounter.size(); i++)
	{
	    int word_index=(pos/per_magnify),offset=(pos%per_magnify);
		int value = (allcounter[i][word_index] >> (offset*magnify_counter_size)) & max_magnify;

		res+=value<<(base_counter_size-1+(i-base_level-1)*(magnify_counter_size-flag_size)+base_level);
		int overflow_flag=((allcounter[i][word_index] >> (offset*magnify_counter_size+magnify_counter_size-flag_size)) &max_flag);
		if(overflow_flag){
            pos=(getpos(pos,i+1)+overflow_flag-1)%(width[i+1]*width[i+1]);
            continue;
		}
		break;
	}
	return res;
}
int PGSS::matrixQuery(int pos){
    int my_word_index=(pos/per_base),counter_offset=(pos%per_base);

    int word_index=(pos/64),offset=(pos%64);
    int value=0;
    for(int j=0;j<base_level;j++){
        value+=(((allcounter[j][word_index] >> offset) & (0x1))<<j);
    }

    value+=(((allcounter[base_level][my_word_index] >> (counter_offset*base_counter_size)) & max_base)<<base_level);

    if ((allcounter[base_level][my_word_index]>>((counter_offset*base_counter_size)+base_counter_size-1))&0x1)
    {
        pos=getpos(pos,base_level+1);
        value+=get_value(pos);
    }
    return value;
}
int PGSS::edgeQuery(string s1, string s2)// s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
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
		for (int j = 0; j<s; j++)
		{

			if ( (value[pos].idx[j] == (index2 | (index1 << 8)))&& (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
			{
				delete []tmp1;
				delete []tmp2;
				return matrixQuery(pos);
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
		delete []tmp1;
		delete []tmp2;
		return 0;
}
void PGSS::nodeSuccessorQuery(string s1, vector<string>&IDs)// query the successors of a node, s1 is the ID of the queried node. results are put in the vector, hash table needed.
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
			for (int j = 0; j < s; ++j)
			{
				if ((value[pos].idx[j]>>8) == i && (value[pos].src[j] == g1))
				{
					     int tmp_g = value[pos].dst[j];
						 int tmp_s = (value[pos].idx[j]&((1<<8)-1));
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
			node = node->next;
			while (node != NULL)
			{
				mapTable.getID(node->key, IDs);
				node=node->next;
			}
		}
		delete []tmp1;
		return;
}
void PGSS::nodePrecursorQuery(string s1, vector<string>&IDs) // same as successor query
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
			for (int j = 0; j < s; ++j)
			{
				if ((value[pos].idx[j]&((1<<8)-1)) == i && (value[pos].dst[j] == g1))
				{
					     int tmp_g = value[pos].src[j];
						 int tmp_s = (value[pos].idx[j]>>8);

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
			node = node->next;
			while (node != NULL)
			{
				if(node->key == k1)
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
//type 0 is for successor query, type 1 is for precusor query
int PGSS::nodeValueQuery(string s1, int type) // s1 is the ID of the queried node, function for node query.
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
				for (int j = 0; j < s; ++j)
				{
					if ((value[pos].idx[j]>>8) == i && (value[pos].src[j] == g1))
					{
						weight += matrixQuery(pos);
					}
				}
			}
			else if (type == 1)/*precursor query*/
			{
				int pos = p1 + k*w;
				for (int j = 0; j < s; ++j)
				{
					if ((value[pos].idx[j]&((1<<8)-1)) == i && (value[pos].dst[j] == g1))
					{
						weight += matrixQuery(pos);
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
			node = node->next;
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
//type 0 is for successor query, type 1 is for precusor query
int PGSS::nodeDegreeQuery(string s1, int type) // s1 is the ID of the queried node, return the in/out degree
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
				for (int j = 0; j < s; ++j)
				{
					if ((value[pos].idx[j]>>8) == i && (value[pos].src[j] == g1))
					{
						 int tmp_g = value[pos].dst[j];
						 int tmp_s = (value[pos].idx[j]&((1<<8)-1));

						 int val=(calHash(tmp_g,tmp_s,k,w)<<f)+tmp_g;


						 degree+=mapTable.countIDnums(val);
					}
				}
			}
			else if (type == 1)/*precursor query*/
			{
				int pos = p1 + k*w;
				for (int j = 0; j < s; ++j)
				{
					if ((value[pos].idx[j]&((1<<8)-1)) == i && (value[pos].dst[j] == g1))
					{
						 int tmp_g = value[pos].src[j];
						 int tmp_s = (value[pos].idx[j]>>8);

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
			node = node->next;
			while (node != NULL)
			{
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
			node = node->next;
			while (node != NULL)
			{
				if (node->key == k1)
				{
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

#endif // PGSS_H_INCLUDED
