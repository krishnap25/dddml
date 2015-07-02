include ../../make/ps_app.mk

CFLAGS  = -O0 -ggdb -Wall -std=c++11 -I./ -I../ $(INCLUDE) $(DMLC_CFLAGS) $(EXTRA_CFLAGS)
LDFLAGS = $(CORE_PATH)/libdmlc.a $(DMLC_LDFLAGS)  $(EXTRA_LDFLAGS)

all: base core bin # ps # tool # test

clean:
	rm -rf bin/*

bin: bin/datasplit

#build/%.o: %.cc
#	@mkdir -p $(@D)
#	$(CXX) $(CFLAGS) -MM -MT build/$*.o $< >build/$*.d
#	$(CXX) $(CFLAGS) -c $< -o $@

bin/%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -MM -MT bin/$*.o $< >bin/$*.d
	$(CXX) $(CFLAGS) -c $< -o $@


#%.pb.cc %.pb.h : %.proto
#	${DEPS_PATH}/bin/protoc --cpp_out=. --proto_path=. $<

bin/datasplit:  ../base/base.a 
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

-include bin/*.d
