#if DISTRIBUTED

#include "dmlc/data.h"
#include "data/row_block.h"

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
	
	/* Step 1: Figure out number of files */
	int nFiles = 1,
		nPartPerFile = 1000,
		nPartToRead = 10,
		mb_size = 1000,
		partID;
	char data_format[] = "libsvm";
	
	real_t probability_of_selecting_one_row = ${SOME_PROB}; //TODO: Calculate.
		
	/* Step 2: Read some of the blocks at random, and sub-sample */
	
	dmlc::data::RowBlockContainer<unsigned> sample;
	
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
				using std::vector;
				using std::shared_ptr;
				using Minibatch = dmlc::data::RowBlockContainer<unsigned>;
				
				auto mb = reader.Value(); //row block
				for (int i = 0; i < mb.size(); ++i)
				{
					//decide whether to add row mb[i] to the sample or not
					
					
				}
			}
		}
	}
	/* Step 3: Localize */
	
	
	/* Step 4: Write to file */


	
}





} //namespace ddmlc

























