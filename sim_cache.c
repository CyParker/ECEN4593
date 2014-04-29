#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <deque>
#include <math.h>
using namespace std;

#define CHUNCKSIZE 16
#define L1_HIT 1
#define L1_MISS 1
#define L2_HIT 5
#define L2_MISS 8
#define L1_L2_TRANSFER 12
#define L2_MAIN_TRANSFER (10+50+(20*(64/CHUNCKSIZE)))

void init_cache(char *configFile);
void LRU_I(int index, int cache_accessed);
void LRU_D(int index, int cache_accessed);
void LRU_2(int index, int cache_accessed);
void L2_Read(uint64_t tag, uint32_t index);
void L2_Write(uint64_t tag, uint32_t index);
void L1I_Read(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index);
void L1D_Write(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index);
void L1D_Read(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index);

struct L1DataConfig{
	int size;
	int assoc;
	int Tshift;
	int Ishift;
	uint32_t  Imask;
};

struct L1DataConfig L1dConfig = {
	8192,
	1,
	13,
	5,
	0x1FE0
};

 struct L1InstConfig{
	int size;
	int assoc;
	int Tshift;
	int Ishift;
	uint32_t  Imask;
};

struct L1InstConfig L1iConfig = {
	8192,
	1,
	13,
	5,
	0x1FE0
};

struct L2CacheConfig{
	int size;
	int assoc;
	int Tshift;
	int Ishift;
	uint32_t  Imask;
};

struct L2CacheConfig L2Config = {
	32768,
	1,
	15,
	6,
	0x7FC0
};

typedef struct CacheBlock{
	int valid;
	int dirty;
	uint64_t tag;
}CACHEBLOCK;


CacheBlock **L1DataCache;
CacheBlock **L1InstCache;
CacheBlock **L2Cache;
std::deque<int> L1dLRUarray[256];
std::deque<int> L1iLRUarray[256];
std::deque<int> L2LRUarray[1024];
unsigned long long execute_time= 0;
unsigned long long read = 0;
unsigned long long write = 0;
unsigned long long inst = 0;

unsigned long long last_time = 0;
unsigned long long read_time = 0;
unsigned long long write_time = 0;
unsigned long long inst_time = 0;
unsigned long long tot_refs = 0;

unsigned long long L1i_req = 0;
unsigned long long L1i_hit = 0;
unsigned long long L1i_kickout = 0;
unsigned long long L1i_dkickout = 0;

unsigned long long L1d_req = 0;
unsigned long long L1d_hit = 0;
unsigned long long L1d_kickout = 0;
unsigned long long L1d_dkickout = 0;

unsigned long long L2_req = 0;
unsigned long long L2_hit = 0;
unsigned long long L2_kickout = 0;
unsigned long long L2_dkickout = 0;

int main(int argc, char *argv[])
{
	char op;
	unsigned long long int total_refs = 0;
	unsigned long long int address;
	uint64_t L1tag;
	uint32_t L1index;
	uint64_t L2tag;
	uint32_t L2index;
	unsigned long long int refTable[10];
	int ref_num = 0;
	unsigned int bytesize;
	uint8_t lowlowbytes;
	int num_bytes;

	
	init_cache((argc > 1) ? argv[1] : NULL);
	
	while(scanf("%c %Lx %d\n", &op, &address,&bytesize) == 3)
	{
		switch(op)
		{
		case 'R':
			read++;
			break;
		case 'W':
			write++;
			break;
		case 'I':
			inst++;
			break;
		}
		//Parse Data
		lowlowbytes = address & (0x3);
		num_bytes = bytesize - (4-lowlowbytes);
		refTable[0] = address & (0xfffffffffffc);
		total_refs++;
		while(num_bytes > 0)
		{
			ref_num++;
			total_refs++;
			refTable[ref_num] = refTable[(ref_num - 1)] + 4;
			num_bytes = num_bytes - 4;
		}
		
		for(;ref_num >= 0; ref_num--)
		{
			//get tag and index
						
			switch(op)
			{
			case 'R':

				//call l1d read function with refTabe[ref_num]
				L1tag = (refTable[ref_num])>> L1dConfig.Tshift;
				L1index = (refTable[ref_num] & L1dConfig.Imask) >> L1dConfig.Ishift;
				
				L2tag = (refTable[ref_num]) >> L2Config.Tshift;
				L2index = (refTable[ref_num] & L2Config.Imask) >> L2Config.Ishift;
				
				L1D_Read(L1tag, L1index, L2tag, L2index);
				read_time += (execute_time - last_time);
				last_time = execute_time;
				break;
			case 'W':

				//call l1d write function with refTabe[ref_num]
				L1tag = (refTable[ref_num]) >> L1dConfig.Tshift;
				L1index = (refTable[ref_num] & L1dConfig.Imask) >> L1dConfig.Ishift;
				
				L2tag = (refTable[ref_num]) >> L2Config.Tshift;
				L2index = (refTable[ref_num] & L2Config.Imask) >> L2Config.Ishift;
				
				L1D_Write(L1tag, L1index, L2tag, L2index);
				write_time += (execute_time - last_time);
				last_time = execute_time;
				break;
			case 'I':

				L1tag = (refTable[ref_num])>> L1iConfig.Tshift;
				L1index = (refTable[ref_num] & L1iConfig.Imask) >> L1iConfig.Ishift;
				
				L2tag = (refTable[ref_num]) >> L2Config.Tshift;
				L2index = (refTable[ref_num] & L2Config.Imask) >> L2Config.Ishift;
				
				L1I_Read(L1tag, L1index, L2tag, L2index);
				inst_time += (execute_time - last_time);
				last_time = execute_time;
			}
		}
		ref_num = 0;
	}
	int numSets;
	numSets = (L1iConfig.size)/(32 * L1iConfig.assoc);
	for(int i = 0; i < numSets; i ++)
	{
		free(L1InstCache[i]);
	}
	free(L1InstCache);
	
	numSets = (L1dConfig.size)/(32 * L1dConfig.assoc);


	for(int i = 0; i < numSets; i ++)
	{
		free(L1DataCache[i]);
	}
	free(L1DataCache);
	
	
	numSets = (L2Config.size)/(64 * L2Config.assoc);

	for(int i = 0; i < numSets; i ++)
	{
		free(L2Cache[i]);
	}	   
	free(L2Cache);
	setvbuf(stdout, NULL, 0, _IONBF);
	printf("-------------------------------------------------------------------------\n");
	printf("              %s         Simulation Results \n", argv[1]);
	printf("-------------------------------------------------------------------------\n\n");
	
	printf("	Memory System:\n");
	printf("		Dcache size = %i : ways = %i : block size = 32\n", L1dConfig.size, L1dConfig.assoc);
	printf("		Icache size = %i : ways = %i : block size = 32\n", L1iConfig.size, L1iConfig.assoc);
	printf("		L2cache size = %i : ways = %i : block size = 32\n", L2Config.size, L2Config.assoc);
	printf("		Memory read time = 50 : chuncksize = %i : chuncktime = 20\n\n", CHUNCKSIZE);
	
	printf("	Execute time = %llu;	Total refs = %llu\n", execute_time, (read+write+inst));
	printf("	Inst refs = %llu;	Data refs = %llu\n\n", inst, (read + write));
	
	float rwi = read + write + inst;
	float p_read = (float) (read)/rwi;
	float p_write = (float) (write)/rwi;
	float p_inst = (float) (inst)/rwi;
	
	printf("	Number of reference types:  [Percentage]\n");
	printf("		Reads  = %llu	[%.2f%%]\n", read, p_read*100);
	printf("		Writes = %llu	[%.2f%%]\n", write,p_write*100);
	printf("		Inst.  = %llu	[%.2f%%]\n", inst, p_inst*100);
	printf("		Total  = %llu\n\n", (read+write+inst));
	
	float pt_read = (float) read_time / (float) execute_time;
	float pt_write = (float) write_time / (float) execute_time;
	float pt_inst = (float) inst_time / (float) execute_time;
	
	printf("	Total cycles for activities:  [Percentage]\n");
	printf("		Reads  = %llu	[%.2f%%]\n", read_time,(float) (100*pt_read));
	printf("		Writes = %llu	[%.2f%%]\n", write_time,(float) (100*pt_write));
	printf("		Inst.  = %llu	[%.2f%%]\n", inst_time,(float) (100*pt_inst));
	printf("		Total  = %llu\n\n", execute_time);
	
	float read_cycles = ((float) (read_time)/ (float) (read));
	float write_cycles = ((float) (write_time)/ (float) (write));
	float inst_cycles = ((float) (inst_time)/ (float) (inst));
	
	float CPI = (float) (read+write+inst+inst)/ (float) (inst);
	
	printf("	Average cycles per activity:\n");
	printf("		Read = %.2f; Write = %.2f; Inst. = %.2f\n",read_cycles, write_cycles,(float) (inst_time/inst));
	printf("	Ideal: Exec. Time = %llu; CPI = %.2f\n", (read+write+inst+inst), CPI);
	CPI = (float) (inst+L1i_req+L1d_req) / (float) (inst);
	printf("	Ideal mis-aligned: Exec. Time = %llu; CPI = %.2f\n\n", (inst+L1i_req+L1d_req), CPI);
	
	float p_hit = (float) (L1i_hit)/(float) (L1i_req);
	float p_miss = 1 - p_hit;
	
	printf("	Memory Level:  L1i\n");
	printf("		Hit Count = %llu	Miss Count = %llu\n", (L1i_hit), L1i_req-L1i_hit);
	printf("		Total Requests = %llu\n", L1i_req);
	printf("		Hit Rate = %.2f%%	Miss Rate = %.2f%%\n",100*p_hit,100*p_miss);
	printf("		Kickouts = %llu; Dirty kickouts = %llu;\n\n", L1i_kickout, L1i_dkickout);
	
	p_hit = (float) (L1d_hit)/(float) (L1d_req);
	p_miss = 1 - p_hit;
	
	printf("	Memory Level:  L1d\n");
	printf("		Hit Count = %llu	Miss Count = %llu\n", (L1d_hit), L1d_req-L1d_hit);
	printf("		Total Requests = %llu\n", L1d_req);
	printf("		Hit Rate = %.2f%%	Miss Rate = %.2f%%\n",100*p_hit,100*p_miss);
	printf("		Kickouts = %llu; Dirty kickouts = %llu;\n\n", L1d_kickout, L1d_dkickout);
	
	p_hit = (float) (L2_hit)/(float) (L2_req);
	p_miss = 1 - p_hit;
	
	printf("	Memory Level:  L2\n");
	printf("		Hit Count = %llu	Miss Count = %llu\n", (L2_hit), L2_req-L2_hit);
	printf("		Total Requests = %llu\n", L2_req);
	printf("		Hit Rate = %.2f%%	Miss Rate = %.2f%%\n",100*p_hit,100*p_miss);
	printf("		Kickouts = %llu; Dirty kickouts = %llu;\n\n", L2_kickout, L2_dkickout);
	
	int acc_cost = log2((L1iConfig.assoc));
	int L1i_cost = (100*((L1iConfig.size)/4096)) +  (100*acc_cost)*((L1iConfig.size)/4096);
	acc_cost = log2((L1dConfig.assoc));
	int L1d_cost = (100*((L1dConfig.size)/4096)) +  (100*acc_cost)*((L1dConfig.size)/4096);
	int L1_cost = L1i_cost + L1d_cost;
	acc_cost = log2((L2Config.assoc));
	int L2_cost = 50 + 50*((L2Config.size)/65536) + 50*acc_cost;
	acc_cost = CHUNCKSIZE/32;
	int Mem_cost = 75 + 100*acc_cost;
	
	printf("	L1 cache cost (Icache $%i) + (Dcache $%i) = %i\n", L1i_cost, L1d_cost, L1_cost); 
	printf("	L2 cache cost = $%i; Memory cost = $%i\n", L2_cost, Mem_cost);
	printf("	Total Cost = %i\n", (L1_cost + Mem_cost + L2_cost));		
	
	return EXIT_SUCCESS;
	
}

/*****************************************************************************/
// \breif 
//
//
//
/*****************************************************************************/
void init_cache(char *configFile)
{
	int numSets;
	FILE *filePtr;
		if(configFile != NULL)
		{
			filePtr = fopen(configFile, "r");
			
			
			//L1 Instruction cache configuration extraction
			fscanf(filePtr, "%i, %i, %i, %i, %x", &(L1iConfig.size), &(L1iConfig.assoc), &(L1iConfig.Tshift),
				   &(L1iConfig.Ishift), &(L1iConfig.Imask));
		}
		
		numSets = (L1iConfig.size)/(32 * L1iConfig.assoc);
		
		//Create Cache Sets
		L1InstCache = (CACHEBLOCK **) malloc(sizeof(CACHEBLOCK *) * numSets);
		
		//Create Cache Ways
		for(int i = 0; i < numSets; i ++)
		{
			L1InstCache[i] = (CACHEBLOCK *) malloc(sizeof(CACHEBLOCK) * L1iConfig.assoc);
		}
		
		if(configFile != NULL)
		{
			//L1 Data cache configuration extraction
			fscanf(filePtr, "%i, %i, %i, %i, %x", &(L1dConfig.size), &(L1dConfig.assoc), &(L1dConfig.Tshift),
				   &(L1dConfig.Ishift), &(L1dConfig.Imask));
		}
		
		numSets = (L1dConfig.size)/(32 * L1dConfig.assoc);
		
		//Create Cache Sets
		L1DataCache = (CACHEBLOCK **) malloc(sizeof(CACHEBLOCK *) * numSets);
		
		//Create Cache Ways
		for(int i = 0; i < numSets; i ++)
		{
			L1DataCache[i] = (CACHEBLOCK *) malloc(sizeof(CACHEBLOCK) * L1dConfig.assoc);
		}
		
		if(configFile != NULL)
		{	   
			//L2 cache configuration extraction
			fscanf(filePtr, "%i, %i, %i, %i, %x,", &(L2Config.size), &(L2Config.assoc), &(L2Config.Tshift),
				   &(L2Config.Ishift), &(L2Config.Imask));
		}
		
		numSets = (L2Config.size)/(64 * L2Config.assoc);
		
		//Create Cache Sets
		L2Cache = (CACHEBLOCK **) malloc(sizeof(CACHEBLOCK *) * numSets);
		
		//Create Cache Ways
		for(int i = 0; i < numSets; i ++)
		{
			L2Cache[i] = (CACHEBLOCK *) malloc(sizeof(CACHEBLOCK) * L2Config.assoc);
		}	   
	 
}
/******************************************************************************/

/******************************************************************************/
void L1D_Read(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index)
{
	int i;
	L1d_req++;
	int way;
	for(i = 0; i < L1dConfig.assoc; i++)
	{
		if(((L1DataCache[index][i]).tag == tag) && ((L1DataCache[index][i]).valid))
		{
			//hit time++
			execute_time += L1_HIT;
			L1d_hit++;
			LRU_D(index, i);
			return;
		}
		if(!((L1DataCache[index][i]).valid))
		{
			break;
		}
	}
	
	execute_time += L1_MISS;
	//L2 read
	
	if(i < L1dConfig.assoc)
	{
		//no kickout
		L2_Read(L2_tag, L2_index);
		L1dLRUarray[index].push_back(i);
		(L1DataCache[index][i]).tag = tag;
		(L1DataCache[index][i]).valid = 1;
		(L1DataCache[index][i]).dirty = 0;
	}
	else
	{
		uint64_t L2W_tag;
		uint32_t L2W_index;
		L1d_kickout++;
		//get LRU way
		way = L1dLRUarray[index].front();
		//if that way is dirty, write it back to L2
		if(((L1DataCache[index][way]).dirty) && ((L1DataCache[index][way]).valid))
		{
			L1d_dkickout++;
			execute_time += L1_L2_TRANSFER;
			//how do we know the L2 tag?
			L2W_tag = ((((L1DataCache[index][way]).tag) << L1dConfig.Tshift) | (index << L1dConfig.Ishift)) >> L2Config.Tshift;
			L2W_index = (((((L1DataCache[index][way]).tag) << L1dConfig.Tshift) | (index << L1dConfig.Ishift)) & L2Config.Imask) >> L2Config.Ishift;
			L2_Write(L2W_tag, L2W_index);
		}
		//update LRU stack
		L2_Read(L2_tag, L2_index);
		LRU_D(index, way);
		
		(L1DataCache[index][way]).tag = tag;
		(L1DataCache[index][way]).valid = 1;
		(L1DataCache[index][way]).dirty = 0;
	}
	
	execute_time += L1_HIT;
}
/******************************************************************************/

/******************************************************************************/
void L1D_Write(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index)
{
	int i;
	L1d_req++;
	int way;
	for(i = 0; i < L1dConfig.assoc; i++)
	{
		if(((L1DataCache[index][i]).tag == tag) && ((L1DataCache[index][i]).valid))
		{
			//hit time++
			execute_time += L1_HIT;
			L1d_hit++;
			LRU_D(index, i);
			L1DataCache[index][i].dirty = 1;
			return;
		}
		if(!((L1DataCache[index][i]).valid))
		{
			break;
		}
	}
	
	execute_time += L1_MISS;

	if(i < L1dConfig.assoc)
	{
		//no kickout
		L2_Read(L2_tag, L2_index);
		L1dLRUarray[index].push_back(i);
		(L1DataCache[index][i]).tag = tag;
		(L1DataCache[index][i]).valid = 1;
		(L1DataCache[index][i]).dirty = 1;
	}
	else
	{
		uint64_t L2W_tag;
		uint32_t L2W_index;
		L1d_kickout++;
		//get LRU way
		way = L1dLRUarray[index].front();
		//if that way is dirty, write it back to L2
		if(((L1DataCache[index][way]).dirty) && ((L1DataCache[index][way]).valid))
		{
			L1d_dkickout++;
			execute_time += L1_L2_TRANSFER;
			//how do we know the L2 tag?
			L2W_tag = ((((L1DataCache[index][way]).tag) << L1dConfig.Tshift) | (index << L1dConfig.Ishift)) >> L2Config.Tshift;
			L2W_index = (((((L1DataCache[index][way]).tag) << L1dConfig.Tshift) | (index << L1dConfig.Ishift)) & L2Config.Imask) >> L2Config.Ishift;
			L2_Write(L2W_tag, L2W_index);
		}
		L2_Read(L2_tag, L2_index);
		//update LRU stack
		LRU_D(index, way);
		
		(L1DataCache[index][way]).tag = tag;
		(L1DataCache[index][way]).valid = 1;
		(L1DataCache[index][way]).dirty = 1;
	}
	
	execute_time += L1_HIT;
	
}
/******************************************************************************/

/******************************************************************************/
void L1I_Read(uint64_t tag, uint32_t index, uint64_t L2_tag, uint32_t L2_index)
{
	int i;
	L1i_req++;
	int way;
	for(i = 0; i < L1iConfig.assoc; i++)
	{
		if(((L1InstCache[index][i]).tag == tag) && ((L1InstCache[index][i]).valid))
		{
			//hit time++
			execute_time += L1_HIT;
			L1i_hit++;
			LRU_I(index, i);

			return;
		}
		if(!((L1InstCache[index][i]).valid))
		{
			break;
		}
	}
	
	execute_time += L1_MISS;
	//L2 read
	
	if(i < L1iConfig.assoc)
	{
		//no kickout
		L2_Read(L2_tag, L2_index);
		L1iLRUarray[index].push_back(i);
		(L1InstCache[index][i]).tag = tag;
		(L1InstCache[index][i]).valid = 1;
		(L1InstCache[index][i]).dirty = 0;
	}
	else
	{
		uint64_t L2W_tag;
		uint64_t L2W_index;
		L1i_kickout++;
		//get LRU way
		way = L1iLRUarray[index].front();
		//if that way is dirty, write it back to L2
		if(((L1InstCache[index][way]).dirty) && ((L1InstCache[index][way]).valid))
		{
			L1i_dkickout++;
			execute_time += L1_L2_TRANSFER;
			//how do we know the L2 tag?
			L2W_tag = ((((L1InstCache[index][way]).tag) << L1iConfig.Tshift) | (index << L1iConfig.Ishift)) >> L2Config.Tshift;
			L2W_index = (((((L1InstCache[index][way]).tag) << L1iConfig.Tshift) | (index << L1iConfig.Ishift)) & L2Config.Imask) >> L2Config.Ishift;
			L2_Write(L2W_tag, L2W_index);
		}
		//update LRU stack
		L2_Read(L2_tag, L2_index);
		LRU_I(index, way);
		(L1InstCache[index][way]).tag = tag;
		(L1InstCache[index][way]).valid = 1;
		(L1InstCache[index][way]).dirty = 0;
	}
	
	execute_time += L1_HIT;
}
/******************************************************************************/

/******************************************************************************/
void L2_Write(uint64_t tag, uint32_t index)
{
	int i;
	L2_req++;
	int way;
	for(i = 0; i < L2Config.assoc; i++)
	{
		if(((L2Cache[index][i]).tag == tag) && ((L2Cache[index][i]).valid))
		{
			execute_time += L2_HIT;
			L2_hit++;
		
			LRU_2(index, i);

			(L2Cache[index][i]).dirty = 1;
			return;
		}
		if(!((L2Cache[index][i]).valid))
		{
			break;
		}
	}
	
	execute_time += (L2_MISS);
	execute_time += (L2_MAIN_TRANSFER);
	if(i < L2Config.assoc)
	{
		//no kickout
		L2LRUarray[index].push_back(i);
		(L2Cache[index][i]).tag = tag;
		(L2Cache[index][i]).valid = 1;
		(L2Cache[index][i]).dirty = 1;
	}
	else
	{
		L2_kickout++;
		//get LRU way
		way = L2LRUarray[index].front();
		//if that way is dirty, write it back to MM
		if(((L2Cache[index][way]).dirty) && ((L2Cache[index][way]).valid))
		{
			L2_dkickout++;
			execute_time += L2_MAIN_TRANSFER;
		}
		//update LRU stack
		LRU_2(index, way);
		
		(L2Cache[index][way]).tag = tag;
		(L2Cache[index][way]).valid = 1;
		(L2Cache[index][way]).dirty = 1;
	}
	execute_time += L2_HIT;
}
/******************************************************************************/

/******************************************************************************/

void L2_Read(uint64_t tag, uint32_t index)
{
	int i;
	L2_req++;
	int way;
	for(i = 0; i < L2Config.assoc; i++)
	{
		if(((L2Cache[index][i]).tag == tag) && ((L2Cache[index][i]).valid))
		{
			execute_time += L2_HIT;
			L2_hit++;
			execute_time += L1_L2_TRANSFER;

			LRU_2(index, i);

			return;
		}
		if(!((L2Cache[index][i]).valid))
		{
			break;
		}
	}
	
	execute_time += L2_MISS;
	execute_time += L2_MAIN_TRANSFER;
	if(i < L2Config.assoc)
	{
		//no kickout
		L2LRUarray[index].push_back(i);
		(L2Cache[index][i]).tag = tag;
		(L2Cache[index][i]).valid = 1;
		(L2Cache[index][i]).dirty = 0;
	}
	else
	{
		L2_kickout++;
		//get LRU way
		way = L2LRUarray[index].front();
		//if that way is dirty, write it back to MM
		if(((L2Cache[index][way]).dirty) && ((L2Cache[index][way]).valid))
		{
			L2_dkickout++;
			execute_time += L2_MAIN_TRANSFER;
		}
		//update LRU stack
		LRU_2(index, way);
		
		(L2Cache[index][way]).tag = tag;
		(L2Cache[index][way]).valid = 1;
		(L2Cache[index][way]).dirty = 0;
	}
	execute_time += L2_HIT;
	execute_time += L1_L2_TRANSFER;
}
/******************************************************************************/

/******************************************************************************/
void LRU_I(int index, int cache_accessed)
{
	std::deque<int> temp;
	int holder;
	
	while(L1iLRUarray[index].back() != cache_accessed)
	{
		temp.push_back(L1iLRUarray[index].back());
		L1iLRUarray[index].pop_back();
	}
	holder = L1iLRUarray[index].back();
	L1iLRUarray[index].pop_back();
	
	while(!temp.empty())
	{
		L1iLRUarray[index].push_back(temp.back());
		temp.pop_back();
	}
	
	L1iLRUarray[index].push_back(holder);
}
/******************************************************************************/

/******************************************************************************/
void LRU_D(int index, int cache_accessed)
{
	std::deque<int> temp;
	int holder;
	
	while(L1dLRUarray[index].back() != cache_accessed)
	{
		temp.push_back(L1dLRUarray[index].back());
		L1dLRUarray[index].pop_back();
	}
	holder = L1dLRUarray[index].back();
	L1dLRUarray[index].pop_back();
	
	while(!temp.empty())
	{
		L1dLRUarray[index].push_back(temp.back());
		temp.pop_back();
	}
	
	L1dLRUarray[index].push_back(holder);
}
/******************************************************************************/

/******************************************************************************/
void LRU_2(int index, int cache_accessed)
{
	std::deque<int> temp;
	int holder;
	
	while(L2LRUarray[index].back() != cache_accessed)
	{
		temp.push_back(L2LRUarray[index].back());
		L2LRUarray[index].pop_back();
	}
	holder = L2LRUarray[index].back();
	L2LRUarray[index].pop_back();
	
	while(!temp.empty())
	{
		L2LRUarray[index].push_back(temp.back());
		temp.pop_back();
	}
	
	L2LRUarray[index].push_back(holder);
}



