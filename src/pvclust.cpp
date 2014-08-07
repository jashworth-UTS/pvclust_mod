#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <cstddef> // for std::ptrdiff_t
#include <algorithm> // for std::stable_sort
#include <stdexcept> // for std::runtime_error
#include <string> // for std::string
#include <list>
#include <vector>
//#include <sstream> // can't use sstream because Rinternal.h macro length() clashes with another
#include <new> // for std::bad_alloc
#include <exception> // for std::exception

// ja for rand()
#include <stdlib.h>
#include <time.h>

typedef std::vector< std::vector<int> > Members;

void get_members(Members & members, const int * const merge){
	size_t const nrow(members.size());
	for(size_t i(0); i<nrow; ++i){
		int ai(merge[i]);
		if(ai<0) {members[i].clear(); members[i].push_back(-1*ai);}
		// note: indices in R's merge table are 1-indexed
		else members[i] = members[ai-1];

		ai = merge[i+nrow];
		if(ai<0) members[i].push_back(-1*ai);
		else members[i].insert(members[i].end(),members[ai-1].begin(),members[ai-1].end());
	}
	for(size_t i(0); i<members.size(); ++i)
		std::sort(members[i].begin(), members[i].end());
}

extern "C" {
  SEXP hc2_split(SEXP merge_) {
    SEXP r = NULL; // return value

    try{

      PROTECT(merge_ = AS_INTEGER(merge_));
      int const n(LENGTH(merge_)), ncol(2);
			int const nrow(n/ncol);
      const int * const merge = INTEGER_POINTER(merge_);
			Members members(nrow);
			get_members(members, merge);
			UNPROTECT(1); // merge_

			SEXP patterns;
			PROTECT(patterns = NEW_CHARACTER(nrow));
			SEXP * p = CHARACTER_POINTER(patterns);

			for(int i(0); i<nrow; ++i){
				std::vector<char> pat(nrow+1,'0');
				for(size_t j(0); j<members[i].size(); ++j) pat[ members[i][j]-1 ] = '1';
				std::string pstr(pat.begin(), pat.end());
				p[i] = mkChar(pstr.c_str());
			}
			PROTECT(r = NEW_LIST(2));
			SET_ELEMENT(r, 0, patterns);

			SEXP memberlist;
			PROTECT(memberlist = NEW_LIST(nrow));

			for(int i(0); i<nrow; ++i){
				size_t msize(members[i].size());
				SEXP member;
				PROTECT(member = NEW_INTEGER(msize));
				int * m = INTEGER_POINTER(member);
				for(size_t j(0); j<msize; ++j) m[j] = members[i][j];
				SET_ELEMENT(memberlist, i, member);
				UNPROTECT(1); // member
			}

			SET_ELEMENT(r, 1, memberlist);
			UNPROTECT(1); // memberlist

      SEXP names;
      PROTECT(names = NEW_CHARACTER(2));
      SET_STRING_ELT(names, 0, COPY_TO_USER_STRING("pattern"));
			SET_STRING_ELT(names, 1, COPY_TO_USER_STRING("member"));
      SET_NAMES(r, names);
			UNPROTECT(1); // names
			UNPROTECT(1); // r
			UNPROTECT(1); // patterns

    } // try
    catch (const std::bad_alloc&) {
      Rf_error( "Memory overflow.");
    }
    catch(const std::exception& e){
      Rf_error( e.what() );
    }
    catch(...){
      Rf_error( "C++ exception (unknown reason)." );
    }

    return r;
  }

  SEXP count_edge_matches(SEXP merge1_, SEXP merge2_){
    SEXP counts = NULL; // return value

//    try{

      PROTECT(merge1_ = AS_INTEGER(merge1_));
      // merge1_ is a two-column matrix of integers from hclust
      int const n1(LENGTH(merge1_)), ncol(2);
			int const nrow1(n1/ncol);
      const int * const merge1 = INTEGER_POINTER(merge1_);
			Members memb1(nrow1);
			get_members(memb1,merge1);
			UNPROTECT(1); // merge1_

      PROTECT(merge2_ = AS_INTEGER(merge2_));
      // merge2_ is a two-column matrix of integers from hclust
      int const n2(LENGTH(merge2_));
			int const nrow2(n2/ncol);
      const int * const merge2 = INTEGER_POINTER(merge2_);
			Members memb2(nrow2);
			get_members(memb2,merge2);
			UNPROTECT(1); // merge2_

			PROTECT(counts = NEW_INTEGER(nrow1));
			int * const countsp = INTEGER_POINTER(counts);

			// speedup: assume unique members and do not re-compare to previously matched ones
			std::vector<bool> matched(nrow2,false);
			for(size_t i1(0); i1<memb1.size(); ++i1){
				size_t const msz1(memb1[i1].size());
				countsp[int(i1)] = 0;
				for(size_t i2(0); i2<memb2.size(); ++i2){
					if(matched[i2]) continue;
					size_t const msz2(memb2[i2].size());
					if(msz1 != msz2) continue;
					bool match(true);
					for(size_t j(0); j<msz2; ++j){
						if(memb1[i1][j] != memb2[i2][j]) {match=false; break;}
					}
					if(match){
						countsp[int(i1)] = 1;
						matched[i2] = true;
					}
				}
			}

			UNPROTECT(1); // counts

//    } // try
//    catch (const std::bad_alloc&) {
//      Rf_error( "Memory overflow.");
//    }
//    catch(const std::exception& e){
//      Rf_error( e.what() );
//    }
//    catch(...){
//      Rf_error( "C++ exception (unknown reason)." );
//    }

    return counts;
  }

  void R_init_pvclust(DllInfo * const info)
  {
    R_CallMethodDef callMethods[]  = {
      {"hc2_split", (DL_FUNC) &hc2_split, 2},
      {"count_edge_matches", (DL_FUNC) &count_edge_matches, 2},
      {NULL, NULL, 0}
    };
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  }

} // extern "C"
