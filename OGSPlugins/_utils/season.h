/*=========================================================================

  Module:    Season

  Implementation of bit.sea Season in C++. 

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef season_h
#define season_h

#include <string>
#include <vector>

#include "TimeObject.h"
#include "TimeInterval.h"

namespace Time {

	static const std::vector<std::string> DEF_START_SEASON{std::string("0101"),std::string("0401"),std::string("0701"),std::string("1001")};
	static const std::vector<std::string> DEF_NAME_SEASON{std::string("winter"),std::string("spring"),std::string("summer"),std::string("fall")};

	class season {

		public:
			inline season();
			inline season(const int reference_year);
			inline season(const std::vector<std::string> &startseason, const std::vector<std::string> nameseason);
			inline season(const int reference_year, std::vector<std::string> &startseason, std::vector<std::string> nameseason);
			inline ~season();

			inline int          get_reference_year() const;
			inline void         set_reference_year(const int reference_year);

			inline void         setseasons(const std::vector<std::string> &startseason, const std::vector<std::string> nameseason);
			inline int          get_seasons_number() const;
			inline TimeInterval get_season_dates(const int season_num, std::string &season_name);
			inline int          findseason(const TimeObject &TO);

		private:
			int numbers_season, reference_year;
			bool alloc;
			TimeObject  *SEASON_LIST;
			std::string *SEASON_LIST_NAME;

			inline void allocate();
			inline void clear();
	};

	inline void season::allocate() {
		if (!this->alloc) {
			this->SEASON_LIST      = new TimeObject[this->numbers_season];
			this->SEASON_LIST_NAME = new std::string[this->numbers_season];
			this->alloc = true;
		}
	}

	inline void season::clear() {
		if (this->alloc) {
			delete [] this->SEASON_LIST;
			delete [] this->SEASON_LIST_NAME;
			this->alloc = false;
		}
	}

	inline season::season() {
		this->numbers_season = 0;
		this->reference_year = 2000;
		this->alloc          = false;
	}

	inline season::season(const int reference_year) {
		this->numbers_season = 0;
		this->reference_year = reference_year;		
		this->alloc          = false;
	}

	inline season::season(const std::vector<std::string> &startseason, const std::vector<std::string> nameseason) {
		this->numbers_season = 0;
		this->reference_year = 2000;
		this->alloc          = false;
		this->setseasons(startseason,nameseason);		
	}

	inline season::season(const int reference_year, std::vector<std::string> &startseason, std::vector<std::string> nameseason) {
		this->numbers_season = 0;
		this->reference_year = reference_year;
		this->alloc          = false;
		this->setseasons(startseason,nameseason);		
	}

	inline season::~season() { this->clear(); }

	inline int season::get_reference_year() const                    { return this->reference_year; }
	inline void season::set_reference_year(const int reference_year) { this->reference_year = reference_year; }

	inline void season::setseasons(const std::vector<std::string> &startseason, const std::vector<std::string> nameseason) {
		/*
        Given two arrays where is defined the date when season start and its name
        the subroutine generate the season list. In input take two arrays of string:

        - Example of startseason : ["1221","0321","0622","0921"]
        - Example of nameseason  : ["winter","spring","summer","fall"]

        In the previous example we shown that winter start to 21 december,
        spring to 21 march, summer to 22 june and fall on 21 september.		
		*/
		this->numbers_season = (int)(startseason.size()); this->allocate();

		if (startseason.size() != nameseason.size()) {
			fprintf(stderr, "ERROR : arrays definitions mistmatch!\n"); exit(-1);
		}

		int ref_year = ( startseason[0] == std::string("0101") ) ? this->reference_year : this->reference_year-1;
		std::string timestr       = std::to_string(ref_year)+startseason[0]+std::string("-00:00:00");
		this->SEASON_LIST[0]      = TimeObject(timestr.c_str(),"%Y%m%d-%H:%M:%S");
		this->SEASON_LIST_NAME[0] = nameseason[0];

		for(int ii=1; ii<this->numbers_season; ++ii) {
			ref_year = this->reference_year;
			std::string timestr        = std::to_string(ref_year)+startseason[ii]+std::string("-00:00:00");;
			this->SEASON_LIST[ii]      = TimeObject(timestr.c_str(),"%Y%m%d-%H:%M:%S");
			this->SEASON_LIST_NAME[ii] = nameseason[ii];
		}
	}

	inline int season::get_seasons_number() const {
		/*
		Return the number of seasons defined in this object
		*/
		return this->numbers_season; 
	}

	inline TimeInterval season::get_season_dates(const int season_num, std::string &season_name) {
		/*
		Given season number, return the range of season dates (start and end)
		and the name of season.		
		*/
		assert(season_num < this->numbers_season);

		TimeObject start_date(this->SEASON_LIST[season_num]);
		TimeObject end_date(this->SEASON_LIST[season_num]);
		
		if (season_num+1 == this->numbers_season) {
			end_date = this->SEASON_LIST[0];
			end_date.increment_year(1);
		} else {
			end_date = this->SEASON_LIST[season_num+1];
		}

		season_name = this->SEASON_LIST_NAME[season_num];

		//TimeInterval TI(start_date,end_date);
		return TimeInterval(start_date,end_date);
	}

	inline int season::findseason(const TimeObject &TO) {
		/*
		Takes a date as input and return the number and name of season where it is in.
		*/
		TimeInterval TI; 
		std::string season_name;
		TimeObject aux(TO);

		int delta_year = this->reference_year - std::stoi(aux.as_string("%Y"));
		
		aux.increment_year( delta_year );
		for (int ii = 0; ii<this->numbers_season; ++ii) {
			TI = this->get_season_dates(ii,season_name);
			if (TI.contains(aux))
				return ii;
		}

		aux.decrement_year( 1 );
		for (int ii = 0; ii<this->numbers_season; ++ii) {
			TI = this->get_season_dates(ii,season_name);
			if (TI.contains(aux))
				return ii;
		}
		return -1;
	}
}

#endif