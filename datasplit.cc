#if DISTRIBUTED

#include "dmlc/data.h"
#include "data/row_block.h"
#include "base/localizer.h"
#include "base/minibatch_iter.h"

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


void ReadFile(const char* featureFile, std::vector<FeaID> *features)
{
	Stream *file = Stream::Create(featureFile, "r", true);
	if (file == NULL)
	{
		std::cerr << "File doesn't exist\n";
		exit(-1);
	}
	else
	{
		file->Read(features);
	}

}

std::vector<FeaID> *Intersect(const std::vector<FeaID> &v1, const std::vector<FeaID> &v2)
{
	std::vector<FeaID> *output = new std::vector<FeaID>();
	for (int i=0,j=0;((i < v1.size()) && (j < v2.size())); )
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
	std::uniform_real_distribution<> dist(0, 1);
	std::uniform_int_distribution<int> dis(0, nPartToRead - 1);
	
	/* Step 1: Figure out number of files */
	int nFiles = 1,
		nPartPerFile = 1000,
		nPartToRead = 10,
		mb_size = 1000,
		partID;
	//char data_format[] = "libsvm";
	
	real_t probability_of_selecting_one_row = (static_cast<real_t> (subsample)) / total_size * nPartPerFile / nPartToRead; //TODO: Check.
	probability_of_selecting_one_row = std::min(1, probability_of_selecting_one_row);
		
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
			dmlc::data::MinibatchIter<FeaID> reader(
				filename, partID, nPartPerFile,
				data_format, mb_size);
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				for (int i = 0; i < mb.size(); ++i)
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
	RowBlock<unsigned> sample1 = sample.GetBlock();
	/* 3.1: read feature file */
	int SomeDefaultStartingValue = 10000; //TODO
	std::vector<FeaID> features;
	features.reserve(SomeDefaultStartingValue);
	ReadFile(featureFile, &features);
	/* 3.2: Get set of features to keep using localizer */
	Localizer <FeaID> lc;
	std::vector<FeaID> *uidx = new std::vector<FeaID>();
	lc.CountUniqIndex<FeaID>(sample1, 4, uidx, NULL);
	/* 3.3: intersect uidx with features */
	std::vector<FeaID> *idx_dict = Intersect(features, *uidx);
	/* 3.4: localize */
	lc.RemapIndex(sample1, *idx_dict, sample_compressed);
	
	
	/* Step 4: Write to file */
	
	//First write idx_dict. Then write compressed sample.

	Stream *output = Stream::Create(outputFile, "w");
	output->Write(*idx_dict);
	sample_compressed->Save(output);
}

} //namespace ddmlc



int main(int argc, char const *argv[])
{
	std::random_device rd;
	std::mt19937_64 rng (rd());	
	return 0;
}





















