#if DISTRIBUTED

#include "dmlc/data.h"
#include "data/row_block.h"
#include "base/localizer.h"

#else

#include "local/data.h"
#include "local/row_block.h"

#endif

#include <cstring>
#include <random>
#include "sample_helper.h"

namespace dddml{
using FeaID = unsigned;
/*
*	Sub-sample data
*	
	TO-DO: 1. get args from conf file
*/


void ReadFile(const char* featureFile, std::vector<FeaID> &features)
{


}

std::vector<FeaID> Intersect(const std::vector<FeaID> &v1, const std::vector<FeaID> &v2)
{

}

void subsample(
	const char* featureFile, //name of feature file
	const char* data_directory,
	const char* outputFile,
	unsigned int subsample_size,
	unsigned int total_size,
	int seed = 0	
)
{
	//RNG: Mersenne-Twister
	std::mt19937_64 rng (seed);	
	std::uniform_real_distribution<> dist(0, 1);
	
	/* Step 1: Figure out number of files */
	int nFiles = 1,
		nPartPerFile = 1000,
		nPartToRead = 10,
		mb_size = 1000,
		partID;
	char data_format[] = "libsvm";
	
	real_t probability_of_selecting_one_row = ${SOME_PROB}; //TODO: Calculate.
		
	/* Step 2: Read some of the blocks at random, and sub-sample */
	
	dmlc::data::RowBlockContainer<FeaID> sample;
	dmlc::data::RowBlockContainer<FeaID> *sample_compressed = new dmlc::data::RowBlockContainer<FeaID>();
	
	for (int fi = 0; fi < nFiles; ++fi)
	{
		std::uniform_int_distribution<int> dis(9, nPartToRead - 1);
		for (int part = 0; part < nPartToRead; ++part)
		{
			partID = dis(rng);
			//TODO: verify filename
			char filename[200];
			std::sprintf(filename, "%s/%d", data_directory, fi);
			dmlc::data::MinibatchIter<FeaID> reader(
				filename, partID, nPartPerFile,
				data_format, mb_size);
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				for (int i = 0; i < mb.size(); ++i)
				{
					//decide whether to add row mb[i] to the sample or not
					if (dist(rng) > 0.5000)
					{
						sample.Push(mb[i]);
					}
					
				}
			}
		}
	}
	/* Step 3: Localize */
	RowBlock<unsigned> sample1 = sample.GetBlock();
	/* 3.1: read feature file */
	int SomeDefaultStartingValue = 10000; //To-Do
	std::vector<FeaID> features(SomeDefaultStartingValue);
	ReadFile(featureFile, features);
	/* 3.2: Get set of features to keep using localizer */
	Localizer <FeaID> lc;
	std::vector<FeaID> *uidx = new std::vector<FeaID>();
	lc.CountUniqIndex<FeaID>(sample1, 4, uidx, NULL);
	/* 3.3: intersect uidx with features */
	std::vector<FeaID> idx_dict = Intersect(features, *uidx);
	/* 3.4: localize */
	lc.RemapIndex(sample1, idx_dict, sample_compressed);
	
	
	/* Step 4: Write to file */
	
	
	
}





} //namespace ddmlc

























