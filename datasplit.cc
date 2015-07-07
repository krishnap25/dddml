#include "base.h"

#if DISTRIBUTED

#include "dmlc/data.h"
#include "data/row_block.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"

#else

#include "local/data.h"
#include "local/row_block.h"

#endif

#include <cstring>
#include <algorithm>
#include <random>
#include "sample_helper.h"

namespace dddml{
using FeaID = unsigned;
using namespace dmlc;
using namespace dmlc::data;

/*
*	Sub-sample data
*	
	TO-DO: 1. get args from conf file
*/


void ReadFile(const char* featureFile, std::vector<FeaID> *features)
{
	dmlc::Stream *file = dmlc::Stream::Create(featureFile, "r", true);
	if (file == NULL)
	{
		std::cerr << "File doesn't exist\n";
		//exit(-1);
	}
	else
	{
		file->Read(features);
	}

}

std::vector<FeaID> *Intersect(std::vector<FeaID> &v1, std::vector<FeaID> &v2)
{
	// features do not exist
	if (v1.size() == 0) return &v2;
	else if (v2.size() == 0) return &v1;
	//else:
	std::vector<FeaID> *output = new std::vector<FeaID>();
	for (unsigned i=0,j=0;((i < v1.size()) && (j < v2.size())); )
	{
		if (v1[i] == v2[j])
		{
			output->push_back(v1[i]);
			++i; ++j;
		}
		else if(v1[i] > v2[j])
		{
			++j;
		}
		else if (v1[i] < v2[j])
		{
			++i;
		}
	}
	return output;

}

void subsample(
	const char* featureFile, //name of feature file
	const char* data_directory,
	const char* outputFile,
	const char *data_format,
	unsigned int subsample_size,
	unsigned int total_size,
	std::mt19937_64 &rng 	
)
{
using real_t = dmlc::real_t;
	
	/* Step 1: Figure out number of files */
	int nFiles = 1,
		nPartPerFile = 1000,
		nPartToRead = 10,
		mb_size = 1000,
		partID;

	std::uniform_real_distribution<> dist(0, 1);
	std::uniform_int_distribution<int> dis(0, nPartToRead - 1);
	
	real_t probability_of_selecting_one_row = (static_cast<real_t> (subsample_size)) / total_size * nPartPerFile / nPartToRead; //TODO: Check.
	probability_of_selecting_one_row = (probability_of_selecting_one_row > 1.0) ? 1.0 : probability_of_selecting_one_row;
		
	/* Step 2: Read some of the blocks at random, and sub-sample */
	
	dmlc::data::RowBlockContainer<FeaID> sample;
	dmlc::data::RowBlockContainer<FeaID> *sample_compressed = new dmlc::data::RowBlockContainer<FeaID>();
	
	for (int fi = 0; fi < nFiles; ++fi)
	{
		for (int part = 0; part < nPartToRead; ++part)
		{
			partID = dis(rng);
			//TODO: verify filename
			char filename[200];
			std::sprintf(filename, "%s/%d", data_directory, fi);
			MinibatchIter<FeaID> reader(
				filename, partID, nPartPerFile,
				data_format, mb_size);
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				for (size_t i = 0; i < mb.size; ++i)
				{
					//decide whether to add row mb[i] to the sample or not
					if (dist(rng) < probability_of_selecting_one_row)
					{
						sample.Push(mb[i]);
					}
					
				}
			}
		}
	}
	/* Step 3: Localize */
	dmlc::RowBlock<unsigned> sample1 = sample.GetBlock();
	/* 3.1: read feature file */
	int SomeDefaultStartingValue = 10000; //TODO
	std::vector<FeaID> features;
	features.reserve(SomeDefaultStartingValue);
	ReadFile(featureFile, &features);
	/* 3.2: Get set of features to keep using localizer */
	dmlc::Localizer <FeaID> lc;
	std::vector<FeaID> *uidx = new std::vector<FeaID>();
	lc.CountUniqIndex<FeaID>(sample1, /*4,*/ uidx, NULL); //include nthreads = 4 for older version
	/* 3.3: intersect uidx with features */
	std::vector<FeaID> *idx_dict = Intersect(features, *uidx);
	/* 3.4: localize */
	lc.RemapIndex(sample1, *idx_dict, sample_compressed);
	
	
	/* Step 4: Write to file */
	
	//First write idx_dict. Then write compressed sample.

	dmlc::Stream *output = dmlc::Stream::Create(outputFile, "w");
	output->Write(*idx_dict);
	sample_compressed->Save(output);
	
}

} //namespace dddml



int main()
{
	using namespace dddml;
	std::random_device rd;
	std::mt19937_64 rng (rd());	
	char featureFile[] = "features.txt";
	char data_directory[] = "./data/mnist.txt";
	char outputFile[] = "./data/mnist.out";
	char data_format[] = "libsvm";
	int subsample_size = 1000;
	int total_size = 60000;

	subsample(featureFile, data_directory, outputFile,data_format, subsample_size, total_size, rng);
	
	return 0;
}





















