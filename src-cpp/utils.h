#include "types.h"

inline bool pos_in_range_sig(int pos1, int pos2, uint32 delta) {
        return abs(pos1 - pos2) <= delta;
}

inline bool pos_in_range(uint32 pos1, uint32 pos2, uint32 delta) {
	return abs(pos1 - pos2) <= delta;
}

inline bool pos_in_range_asym(uint32 pos1, uint32 pos2, uint32 delta1, uint32 delta2) {
	uint32 delta11 = pos2 > delta1 ? delta1 : pos2;
	return pos1 > pos2 - delta11 && pos1 < pos2 + delta2;
}
