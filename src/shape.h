#ifndef SHAPE_H
#define SHAPE_H

#include <cinttypes>

class shape
{
protected:
	uint8_t tType = 0; // defines tissue type (linked to optical properties)
	uint8_t priority = 0; // defines building priority (ie when we overwrite other volumes)
public:
	void set_priority(const uint8_t _priority);
	uint8_t get_priority() const {return priority;};

	void set_tType(const uint8_t _tType);
	uint8_t get_tType() const {return tType;};
};

#endif