/*=========================================================================

  Module:    Time List

  Implementation of bit.sea TimeList in C++. 

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef TimeList_h
#define TimeList_h

#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cmath>

#include "TimeObject.h"
#include "TimeInterval.h"
#include "TimeRequestors.h"

#define TIME_OK   0
#define TIME_ERR -1

namespace Time {

	typedef std::vector<Requestor *> REQ_LIST;
	typedef std::vector< std::pair<TimeObject,std::vector<int>> > COUPLED_LIST;

	inline TimeInterval computeTimeWindow(std::string &freqString, TimeObject &currentDate) {
		int year, month, day;
		TimeInterval TI;
		if (freqString == std::string("daily")) {
			day   = std::stoi(currentDate.as_string("%d"));
			month = std::stoi(currentDate.as_string("%m"));
			year  = std::stoi(currentDate.as_string("%Y"));
			TI    = Daily_req(year,month,day).interval();
		}
		if (freqString == std::string("weekly")) {
			day   = std::stoi(currentDate.as_string("%d"));
			month = std::stoi(currentDate.as_string("%m"));
			year  = std::stoi(currentDate.as_string("%Y"));
			TI    = Weekly_req(year,month,day).interval();
		}
		if (freqString == std::string("monthly")) {
			month = std::stoi(currentDate.as_string("%m"));
			year  = std::stoi(currentDate.as_string("%Y"));
			TI    = Monthly_req(year,month).interval();
		}
		if (freqString == std::string("yearly")) {
			year  = std::stoi(currentDate.as_string("%Y"));
			TI    = Yearly_req(year).interval();
		}
		if (freqString.substr(0,5) == std::string("days=")) {
			day   = std::stoi(currentDate.as_string("%d"));
			month = std::stoi(currentDate.as_string("%m"));
			year  = std::stoi(currentDate.as_string("%Y"));
			int ndays = std::stoi(freqString.substr(5));
			TI    = Interval_req(year,month,day,12,ndays).interval();
		}
		return TI;
	}

	inline void deallocList(REQ_LIST req_list) {
		for (Requestor *req : req_list)
			delete req;
		req_list.resize(0);
	}
	
	class TimeList {
		public:
			inline TimeList();
			inline TimeList(const TimeList &TL);
			inline TimeList(const TimeObjectList &datelist);
			inline TimeList(const TimeObjectList &datelist, const char *forceFrequency);
			inline TimeList(const TimeObject &start, const TimeObject &end, const char *delta);
			inline TimeList(const TimeObject &start, const TimeObject &end, const char *delta, const char *forceFrequency);
			// Implementation not needed
//			inline TimeList(const TimeInterval *TL, const char *inputdir, const char *searchstring, 
//				const char *filtervar="None", const char *prefix="ave.", const char *dateformat="%Y%m%d-%H:%M:%S",
//				const int hour=12, const char *forceFrequency="None");

			inline int            len() const       { return nTimes; }
			inline TimeObjectList list() const      { return TOL; }
			inline TimeInterval   interval() const  { return TI; }
			inline std::string    frequency() const { return inputFrequency; }

			inline int      select(Requestor *req);
			inline int      select(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);

			inline REQ_LIST getDailyList();
			inline REQ_LIST getWeeklyList(const int weekday);
			inline REQ_LIST getMonthlist(const bool extrap);
			inline REQ_LIST	getSeasonList(season &s_obj, const bool extrap);
			inline REQ_LIST getYearlist();
			inline REQ_LIST getDecadalList();
			inline REQ_LIST getOwnList();
			inline REQ_LIST getSpecificIntervalList(const char *startstr, const char *fmt, const int days);

			inline int          find(const TimeObject &TO)         { return TOL.find(TO); }
			inline void         merge(TimeObjectList &datelist);
			inline void         merge(TimeList &TL);
			inline COUPLED_LIST couple(TimeObjectList &datetimelist);

			inline TimeObject  operator[](int i) const             { return (i>=0) ? TOL[i] : TOL[nTimes+i]; }
			inline TimeObject &operator[](int i)                   { return (i>=0) ? TOL[i] : TOL[nTimes+i]; }
			inline TimeList   &operator=(const TimeList &TL)       { nTimes = TL.nTimes; TOL = TL.TOL; TI = TL.TI; inputFrequency = TL.inputFrequency; return (*this); }

			class iterator {
				public:
					typedef TimeObject              value_type;
					typedef std::ptrdiff_t          difference_type;
					typedef TimeObject*             pointer;
					typedef TimeObject&             reference;
					typedef std::input_iterator_tag iterator_category;

					inline iterator() : TL(nullptr), i(0) {}
					inline iterator(TimeList* TLi, int ii) : TL(TLi), i(ii) {}

					inline       TimeObject& operator*()            { return (*TL)[i]; }
					inline const TimeObject& operator*() const      { return (*TL)[i]; }

					inline       iterator& operator++()             { ++i; return *this; }
					inline       iterator& operator--()             { --i; return *this; }
					inline       iterator  operator++(int)          { iterator r(*this); ++i; return r; }
					inline       iterator  operator--(int)          { iterator r(*this); --i; return r; }

					inline       iterator& operator+=(int n)        { i += n; return *this; }
					inline       iterator& operator-=(int n)        { i -= n; return *this; }

					inline       iterator  operator+(int n) const   { iterator r(*this); return r += n; }
					inline       iterator  operator-(int n) const   { iterator r(*this); return r -= n; }

					inline difference_type operator-(iterator const& r) const { return i - r.i; }

					inline bool operator<(iterator const& r)  const  { return i <  r.i; }
					inline bool operator<=(iterator const& r) const  { return i <= r.i; }
					inline bool operator>(iterator const& r)  const  { return i >  r.i; }
					inline bool operator>=(iterator const& r) const  { return i >= r.i; }
					inline bool operator!=(const iterator &r) const  { return i != r.i; }
					inline bool operator==(const iterator &r) const  { return i == r.i; }

					inline int  ind()                                { return i; }

				private:
					TimeList* TL;
					int  i;
			};

			inline iterator begin()                               { return iterator{this,0}; }
			inline iterator end()                                 { return iterator{this,nTimes}; }

		private:
			int nTimes;
			std::string inputFrequency;
			TimeObjectList TOL;
			TimeInterval TI;

			inline std::string searchFrequency();
			inline int generaltimeselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generalhourselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generaldayselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generalweekselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generalmonthselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generalclimseasonselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
			inline int generalclimmonthselector(Requestor *req, std::vector<int> &selection, std::vector<double> &weights);
	};

	inline TimeList::TimeList() { 
		this->nTimes = 0; 
	}

	inline TimeList::TimeList(const TimeList &TL) {
		this->nTimes         = TL.nTimes;
		this->TOL            = TL.TOL;
		this->TI             = TL.TI;
		this->inputFrequency = TL.inputFrequency;
	}

	inline TimeList::TimeList(const TimeObjectList &datelist) {
		this->nTimes         = datelist.len();
		this->TOL            = TimeObjectList(datelist);
		this->TI             = TimeInterval(datelist[0],datelist[-1]);
		if (this->nTimes > 1) {
			this->inputFrequency = searchFrequency();
		}
	}

	inline TimeList::TimeList(const TimeObjectList &datelist, const char *forceFrequency) {
		this->nTimes         = datelist.len();
		this->TOL            = TimeObjectList(datelist);
		this->TI             = TimeInterval(datelist[0],datelist[-1]);
		this->inputFrequency = std::string(forceFrequency);
	}

	inline TimeList::TimeList(const TimeObject &start, const TimeObject &end, const char *delta) {
		this->TOL            = TimeObjectList(start,end,delta);
		this->nTimes         = this->TOL.len();
		this->TI             = TimeInterval(this->TOL[0],this->TOL[-1]);
		if (this->nTimes > 1) {
			this->inputFrequency = searchFrequency();
		}
	}

	inline TimeList::TimeList(const TimeObject &start, const TimeObject &end, const char *delta, const char *forceFrequency) {
		this->TOL            = TimeObjectList(start,end,delta);
		this->nTimes         = this->TOL.len();
		this->TI             = TimeInterval(this->TOL[0],this->TOL[-1]);
		this->inputFrequency = std::string(forceFrequency);
	}

	inline int TimeList::select(Requestor *req) {
		/*
			Used to select a single time (or file) regardless to time aggregation

			index = select_one(requestor)
			Returned values:
				- an integer index indicating to access selected times (or files)

			Doesn't work with: Clim_day, Clim_month and Clim_season
		*/

		for(int ii=0;ii<this->nTimes;++ii) {
			if (req->contains(this->TOL[ii]))
				return ii;
		}
		return TIME_ERR;
	}

	inline int TimeList::select(Requestor *req, std::vector<int> &selection, std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		if (req->type() == std::string("hourly")  || req->type() == std::string("clim_hour"))
			return this->generalhourselector(req,selection,weights);
		if (req->type() == std::string("daily")   || req->type() == std::string("clim_day"))
			return this->generaldayselector(req,selection,weights);
		if (req->type() == std::string("weekly"))
			return this->generalweekselector(req,selection,weights);
		if (req->type() == std::string("monthly"))
			return this->generalmonthselector(req,selection,weights);
		if (req->type() == std::string("clim_season"))
			return this->generalclimseasonselector(req,selection,weights);
		if (req->type() == std::string("clim_month"))
			return this->generalclimmonthselector(req,selection,weights);

		return this->generaltimeselector(req,selection,weights);
	}

	inline REQ_LIST TimeList::getDailyList() {
		/*
			Tested only for mooring case, interval = 3 hours
		*/
		REQ_LIST req_list;
		TimeObject starting_centered_day = TimeObject(this->TOL[0].as_string("%Y%m%d"),"%Y%m%d");
		TimeObjectList TL(starting_centered_day,this->TOL[-1],"days=1");

		std::vector<int> indexes;
		std::vector<double> weights;

		for(int ii=0; ii<TL.len(); ++ii) {
			int day   = std::stoi(TL[ii].as_string("%d"));
			int month = std::stoi(TL[ii].as_string("%m"));
			int year  = std::stoi(TL[ii].as_string("%Y"));
			Daily_req D(year,month,day);

			this->generaldayselector(&D,indexes,weights);
			if (indexes.size() > 0)
				req_list.push_back( new Daily_req(year,month,day) );
		}

		return req_list;
	}

	inline REQ_LIST TimeList::getWeeklyList(const int weekday) {
		/*
			WeekList = getWeekList(weekday)
			Weekday is the same of the datetime.isoweekay() method.
			Monday == 1 ... Sunday == 7

			Returns an ordered list of requestors, Weekly_req objects.
		*/
		int PossibleShifts[] = {3,2,1,0,-1,-2,-3}; int index = 0;
		TimeObject starting_centered_day = TimeObject( this->TI.get_start_time() );
		for(int day=weekday-3, count=0; day<weekday+4; ++day,++count) {
			int interested_weekday = (int)(day%7);
			if (interested_weekday == 0) interested_weekday = 7;
			if (interested_weekday == starting_centered_day.isoweekday()) index = count;
		}

		starting_centered_day.increment_day(PossibleShifts[index]);
		TimeObjectList TL(starting_centered_day,this->TOL[-1],"days=7");

		REQ_LIST req_list;
		std::vector<int> indexes;
		std::vector<double> weights;

		for(int ii=0; ii<TL.len(); ++ii) {
			int day   = std::stoi(TL[ii].as_string("%d"));
			int month = std::stoi(TL[ii].as_string("%m"));
			int year  = std::stoi(TL[ii].as_string("%Y"));
			Weekly_req W(year,month,day);

			this->generalweekselector(&W,indexes,weights);
			if (indexes.size() > 0)
				req_list.push_back( new Weekly_req(year,month,day) );
		}

		return req_list;
	}

	inline REQ_LIST TimeList::getMonthlist(const bool extrap=false) {
		/*
			Returns an ordered list of requestors, Monthly_req objects.
			By setting extrap=True, this method extrapolates out of the indicated period (Starttime,EndTime)
			Example: if the input is weekly, centered in 20120301, and Starttime=20120301,
			then extrapolation will return also 201202, because of the part of the week before the centered time.
		*/
		std::vector<std::pair<int,int>> month_list;
		month_list.push_back( std::make_pair(std::stoi(this->TOL[0].as_string("%Y")),std::stoi(this->TOL[0].as_string("%m"))) );

		for (int ii=1; ii<this->nTimes; ++ii) {
			std::pair<int,int> newmonth(std::stoi(this->TOL[ii].as_string("%Y")),std::stoi(this->TOL[ii].as_string("%m")));
			if (std::find(month_list.begin(),month_list.end(),newmonth) == month_list.end())
				month_list.push_back( newmonth );
		}

		if (!extrap) {
			std::vector<std::pair<int,int>> month_list_red;
			TimeObject firstMonth(this->TI.get_start_time().as_string("%Y%m01-00:00:00"),"%Y%m%d-%H:%M:%S");
			TimeObject  lastMonth(this->TI.get_end_time().as_string("%Y%m01-00:00:00"),"%Y%m%d-%H:%M:%S");
			for(int ii=0; ii<month_list.size(); ++ii) {
				char buff[256];
				std::sprintf(buff,"%d%02d01-00:00:00",month_list[ii].first,month_list[ii].second);
				TimeObject firstOfMonth(buff,"%Y%m%d-%H:%M:%S");
				if (firstOfMonth >= firstMonth && firstOfMonth <= lastMonth)
					month_list_red.push_back( month_list[ii] );
			}
			month_list = month_list_red;
		}

		REQ_LIST req_list;
		for(int ii=0; ii<month_list.size(); ++ii)
			req_list.push_back( new Monthly_req(month_list[ii].first,month_list[ii].second) );

		return req_list;
	}

	inline REQ_LIST TimeList::getSeasonList(season &s_obj, const bool extrap=false) {
		/*
			Returns an ordered list of requestors, Season_req objects.
			By setting extrap=True, this method extrapolates out of the indicated period (Starttime,EndTime)
			Example: if the input is weekly, centered in 20120301, and Starttime=20120301,
			then extrapolation will return also 201202, because of the part of the week before the centered time.
		*/
		std::vector<std::pair<int,int>> season_list;
		season_list.push_back( std::make_pair(std::stoi(this->TOL[0].as_string("%Y")),s_obj.findseason(this->TOL[0])) );

		for (int ii=1; ii<this->nTimes; ++ii) {
			std::pair<int,int> newseason( std::stoi(this->TOL[ii].as_string("%Y")),s_obj.findseason(this->TOL[ii]) );
			if ( std::find(season_list.begin(),season_list.end(),newseason) == season_list.end() )
				season_list.push_back( newseason );
		}

		if (!extrap) {
			std::vector<std::pair<int,int>> season_list_red;
			TimeObject first = this->TI.get_start_time();
			TimeObject last  = this->TI.get_end_time();
			Season_req firstSeason(std::stoi(first.as_string("%Y")),s_obj.findseason(first),s_obj);
			Season_req lastSeason(std::stoi(last.as_string("%Y")),s_obj.findseason(last),s_obj);
			for(int ii=0; ii<season_list.size(); ++ii) {
				Season_req req(season_list[ii].first,season_list[ii].second,s_obj);
				first  = firstSeason.interval().get_start_time();
				last   = lastSeason.interval().get_end_time();
				TimeObject ffirst = req.interval().get_start_time();
				TimeObject llast  = req.interval().get_end_time();
				if ( ffirst >= first && llast <= last)
					season_list_red.push_back( season_list[ii] );
			}
			season_list = season_list_red;
		}

		REQ_LIST req_list;
		for(int ii=0; ii<season_list.size(); ++ii)
			req_list.push_back( new Season_req(season_list[ii].first,season_list[ii].second,s_obj) );

		return req_list;
	}

	inline REQ_LIST TimeList::getYearlist() {

		std::vector<int> year_list;
		for (int ii=0; ii<this->nTimes; ++ii) {
			int year = std::stoi( this->TOL[ii].as_string("%Y") );
			if (std::find(year_list.begin(),year_list.end(),year) == year_list.end())
				year_list.push_back(year);
		}

		REQ_LIST req_list;
		for (int ii=0; ii<year_list.size();++ii) {
			req_list.push_back( new Yearly_req(year_list[ii]) );
		}

		return req_list;
	}

	inline REQ_LIST TimeList::getDecadalList() {
		REQ_LIST req_list;
		std::fprintf(stderr, "Not implemented!\n");
		return req_list;
	}

	inline REQ_LIST TimeList::getOwnList() {
		/*
			Not useful for time aggregation, but to get requestors in order to match with observations.

			Must call deallocList to clean memory after the vector is used.
		*/
		REQ_LIST req_list;

		if (this->inputFrequency == std::string("daily")) {
			for(int ii = 0; ii <this->nTimes; ++ii) {
				int year  = std::stoi( this->TOL[ii].as_string("%Y") );
				int month = std::stoi( this->TOL[ii].as_string("%m") );
				int day   = std::stoi( this->TOL[ii].as_string("%d") );
				req_list.push_back( new Daily_req(year,month,day) );
			}
		}
		if (this->inputFrequency == std::string("10days")) {
			for(int ii = 0; ii <this->nTimes; ++ii) {
				int year  = std::stoi( this->TOL[ii].as_string("%Y") );
				int month = std::stoi( this->TOL[ii].as_string("%m") );
				int day   = std::stoi( this->TOL[ii].as_string("%d") );
				req_list.push_back( new Interval_req(year,month,day,12,10) );
			}			
		}
		if (this->inputFrequency == std::string("weekly")) {
			for(int ii = 0; ii <this->nTimes; ++ii) {
				int year  = std::stoi( this->TOL[ii].as_string("%Y") );
				int month = std::stoi( this->TOL[ii].as_string("%m") );
				int day   = std::stoi( this->TOL[ii].as_string("%d") );
				req_list.push_back( new Weekly_req(year,month,day) );
			}			
		}
		if (this->inputFrequency == std::string("monthly")) {
			for(int ii = 0; ii <this->nTimes; ++ii) {
				int year  = std::stoi( this->TOL[ii].as_string("%Y") );
				int month = std::stoi( this->TOL[ii].as_string("%m") );
				req_list.push_back( new Monthly_req(year,month) );
			}			
		}
		if (this->inputFrequency == std::string("yearly")) {
			for(int ii = 0; ii <this->nTimes; ++ii) {
				int year  = std::stoi( this->TOL[ii].as_string("%Y") );
				req_list.push_back( new Yearly_req(year) );
			}			
		}

		if (req_list.size() == 0)
			std::fprintf(stderr, "Not implemented!\n");

		return req_list;
	}

	inline REQ_LIST TimeList::getSpecificIntervalList(const char *startstr, const char *fmt, const int days = 10) {
		/*
			Useful in case of 10 days average, for example
		*/
		TimeObject starttime(startstr,fmt);
		TimeObject endtime = this->TI.get_end_time();
		char buff[256];
		std::sprintf(buff,"days=%d",days);
		TimeObjectList dl(starttime,endtime,buff);

		REQ_LIST req_list;
		for (TimeObject dateobj : dl) {
			int year  = std::stoi(dateobj.as_string("%Y"));
			int month = std::stoi(dateobj.as_string("%m"));
			int day   = std::stoi(dateobj.as_string("%d"));
			int hour  = std::stoi(dateobj.as_string("%H"));
			req_list.push_back( new Interval_req(year,month,day,hour,days) );
		}

		return req_list;
	}

	inline std::string TimeList::searchFrequency() {
		/*
			Returns strings: "daily", "weekly", "monthly", "hourly", "10days"
		*/
		if (this->TOL.len() < 2) {
			std::printf("Frequency cannot be calculated in between %s\n",this->TI.as_string("%Y%m%d").c_str());
			return std::string("None");
		}

		// Compute time difference in seconds
		std::vector<time_t> diffs(this->nTimes-1);
		for (int ii=0; ii<this->nTimes-1; ++ii)
			diffs[ii] = this->TOL[ii+1] - this->TOL[ii];

		// Sort the differences
		std::sort(diffs.begin(),diffs.end());

		// Count the moda
		time_t number = diffs[0], moda = diffs[0];
		int count = 1, count_moda = 1;

		for (int ii=1; ii<this->nTimes-1; ++ii) {
			if (diffs[ii] == number)  {
				// Count occurences of the same number
				++count; 
			} else {
				// This is a new number
				if (count > count_moda) {
					count_moda = count;
					moda       = number;
				}
				count  = 1;
				number = diffs[ii]; 
			}
		}

		float days = (float)(moda/3600./24.);

		// Output
		if (days == 1.)                          return std::string("daily");
		if (days > 6.   && days < 8.)            return std::string("weekly");
		if (days > 26.  && days < 32.)           return std::string("monthly");
		if (days < 1.)                           return std::string("hourly");
		if (days > 364. && days < 367.)          return std::string("yearly");
		if (std::fabs((int)(days) - days) < 0.1) return std::string("days=")+std::to_string((int)(days)); 
		if (days == 10.)                         return std::string("10days");
		if (days > 1.   && days < 7.)            return std::string("None");

		std::fprintf(stderr, "Oopsie! this shouldn't have happened!\n");
		exit(-1);
	}

	inline int TimeList::generaltimeselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		// Separate by input frequency
		if (this->inputFrequency == std::string("daily")) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				if ( req->contains(this->TOL[ii]) ) {
					selection.push_back(ii);
					weights.push_back(1.);
				}
			}
			return TIME_OK;
		}
		if (this->inputFrequency == std::string("weekly")  ||
			this->inputFrequency == std::string("monthly") ||
			this->inputFrequency == std::string("yearly")  ||
			this->inputFrequency == std::string("10days") ) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				TimeInterval T1 = computeTimeWindow(this->inputFrequency,this->TOL[ii]);
				TimeInterval T2 = TimeInterval(req->interval());
				time_t weight = T1.overlapTime(T2);
				if (weight > 0) {
					selection.push_back(ii);
					weights.push_back((double)(weight));
				}
			}
			return TIME_OK;
		}
		return TIME_ERR;
	}

	inline int TimeList::generalhourselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		if ( this->inputFrequency != std::string("hourly") ) 
			return TIME_ERR;

		for (int ii=0; ii<this->nTimes; ++ii) {
			if (req->contains(this->TOL[ii])) {
				selection.push_back(ii);
				weights.push_back(1.);
			}
		}

		return TIME_OK;
	}

	inline int TimeList::generaldayselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		// hourly values are treated as instantaneous values, not time averages
		if ( !(this->inputFrequency == std::string("hourly") ||
		       this->inputFrequency == std::string("daily")) ) // it does not matter how many hours
			return TIME_ERR;

		for (int ii=0; ii<this->nTimes; ++ii) {
			if (req->contains(this->TOL[ii])) {
				selection.push_back(ii);
				weights.push_back(1.);
			}
		}

		return TIME_OK;
	}

	inline int TimeList::generalweekselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		// Only works with daily data
		if ( this->inputFrequency != std::string("daily") )
			return TIME_ERR;

		for (int ii=0; ii<this->nTimes; ++ii) {
			if (req->contains(this->TOL[ii])) {
				selection.push_back(ii);
				weights.push_back(1.);
			}
		}

		return TIME_OK;
	}

	inline int TimeList::generalmonthselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		if (this->inputFrequency == std::string("daily") ||
			this->inputFrequency == std::string("hourly")) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				if (req->contains(this->TOL[ii])) {
					selection.push_back(ii);
					weights.push_back(1.);
				}
			}
			return TIME_OK;
		}
		if (this->inputFrequency == std::string("weekly")) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				TimeInterval T1 = computeTimeWindow(this->inputFrequency,this->TOL[ii]);
				TimeInterval T2 = TimeInterval( req->interval() );
				time_t weight = T1.overlapTime(T2);
				if (weight > 0) {
					selection.push_back( ii );
					weights.push_back( (double)(weight) );
				}
			}
			return TIME_OK;
		}
		if (this->inputFrequency.substr(0,5) == std::string("days=")) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				TimeInterval T1 = computeTimeWindow(this->inputFrequency,this->TOL[ii]);
				TimeInterval T2 = TimeInterval( req->interval() );
				time_t weight = T1.overlapTime(T2);
				if (weight > 0) {
					selection.push_back( ii );
					weights.push_back( (double)(weight) );
				}
			}
			return TIME_OK;
		}
		if (this->inputFrequency == std::string("monthly")) {
			// No time aggregation
			for (int ii=0; ii<this->nTimes; ++ii) {
				if (req->contains(this->TOL[ii])) {
					selection.push_back(ii);
					weights.push_back(1.);
				}
			}
			return TIME_OK;
		}
		return TIME_ERR;
	}

	inline int TimeList::generalclimseasonselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {

		// Separate by input frequency
		if (this->inputFrequency == std::string("daily")) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				if ( req->contains(this->TOL[ii]) ) {
					selection.push_back(ii);
					weights.push_back(1.);
				}
			}
			return TIME_OK;
		}
		if (this->inputFrequency == std::string("weekly")  ||
			this->inputFrequency == std::string("monthly") ||
			this->inputFrequency == std::string("yearly")  ||
			this->inputFrequency == std::string("10days") ) {
			for (int ii=0; ii<this->nTimes; ++ii) {
				TimeInterval T1 = computeTimeWindow(this->inputFrequency,this->TOL[ii]);
				TimeInterval T2 = TimeInterval(req->interval());
				int ref_year = std::stoi(req->interval().get_start_time().as_string("%Y"));
				int year     = std::stoi(this->TOL[ii].as_string("%Y"));
				int d_years  = year - ref_year;
				TimeObject starttime(T2.get_start_time()); starttime.increment_year(d_years);
				TimeObject endtime(T2.get_end_time());     endtime.increment_year(d_years);
				T2 = TimeInterval(starttime,endtime);
				time_t weight = T1.overlapTime(T2);
				if (weight > 0) {
					selection.push_back(ii);
					weights.push_back((double)(weight));
				}
			}
			return TIME_OK;
		}
		return TIME_ERR;
	}

	inline int TimeList::generalclimmonthselector(Requestor *req, std::vector<int> &selection, 
		std::vector<double> &weights) {
		/*
			Method for time aggregation
			indexes, weights = select(requestor)
			Returned values:
				- a list of indexes (integers) indicating to access selected times (or files)
				- an array of weights
		*/
		int month = std::stoi( req->interval().get_start_time().as_string("%m") );

		REQ_LIST yearlist = this->getYearlist();
		for(int ii=0; ii<yearlist.size();++ii) {
			int year = std::stoi( yearlist[ii]->interval().get_start_time().as_string("%Y") );
			Monthly_req req(year,month);
			if (this->select(&req,selection,weights) == TIME_ERR) {
				deallocList(yearlist);
				return TIME_ERR;
			}
		}
		deallocList(yearlist);
		return TIME_OK;
	}

	inline COUPLED_LIST TimeList::couple(TimeObjectList &datetimelist) {
		/*
			Performs the association between two timelists:
				- a regular one, usual result of model
				- an irregular one, as provided from measurements

			Input:
				* datetimelist * a list of datetime objects, usually the irregular times of profiles objects

			Returns:
				* Coupled_List * a list of tuples (datetime object, list_of_indices)
			where
				- the datetime object is an element of self.Timelist (model regular times) concerned by datetimelist
				- list_of_indices is a list of integers such that :
				datetimelist[list_of_indeces] belong to the datetime object corresponding timeinterval
		*/
		COUPLED_LIST coupled_list;

		REQ_LIST req_list = this->getOwnList();
		for (int ii=0; ii<req_list.size(); ++ii) {
			std::vector<int> list_of_ind;
			for (int jj=0; jj<datetimelist.len(); ++jj) {
				if (req_list[ii]->contains(datetimelist[ii]))
					list_of_ind.push_back(jj);
			}
			if(list_of_ind.size() > 0)
				coupled_list.push_back( std::make_pair(this->TOL[ii],list_of_ind) );
		}

		deallocList(req_list);
		return coupled_list;
	}

	inline void TimeList::merge(TimeObjectList &datelist) {
		this->TOL            = this->TOL.merge(datelist);
		this->nTimes         = this->TOL.len();
		this->TI             = TimeInterval(this->TOL[0],this->TOL[-1]);
		if (this->nTimes > 1) {
			this->inputFrequency = searchFrequency();
		}
	}

	inline void TimeList::merge(TimeList &TL) { 
		TimeList aux(TL.list());
		if (aux.frequency() == TL.frequency()) 
			this->merge(aux); 
	}
}

#endif