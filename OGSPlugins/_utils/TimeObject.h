/*=========================================================================

  Module:    Time Object

  Class to manage struct tm. 

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef TimeObject_h
#define TimeObject_h

#include <ctime>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

namespace Time {
	
	class TimeObject {

		public:
			inline TimeObject()                                                            { now(); }
			inline TimeObject(tm &in)                                                      { set_tm(in); }
			inline TimeObject(const TimeObject &TO)                                        { set_tm(TO.ltm); }
			inline TimeObject(const char *timestr, const char *fmt)                        { now(); from_string(timestr,fmt); }
			inline TimeObject(const std::string &timestr, const char *fmt)                 { now(); from_string(timestr.c_str(),fmt); }
			inline TimeObject(const std::string &timestr, const std::string &fmt)          { now(); from_string(timestr,fmt); }

			inline tm          get_tm() const                                              { return ltm; }
			inline void        set_tm(const tm &in)                                        { ltm = in; update(); }
			inline time_t      epoch() const                                               { return t; }
			inline int         isoweekday()                                                { return std::stoi(as_string("%u")); }

			inline void        from_string(const char *str, const char *fmt)               { strptime(str,fmt,&ltm); update(); }
			inline void        from_string(const std::string &str, const std::string &fmt) { from_string(str.c_str(),fmt.c_str()); }
			inline std::string as_string(const char *fmt)                                  { char buff[256]; std::strftime(buff,256,fmt,&ltm); return std::string(buff); }
			inline std::string as_string(std::string &fmt)                                 { return as_string(fmt.c_str()); }

			inline void        increment_sec(const int i)                                  { ltm.tm_sec  += i; update(); }
			inline void        increment_min(const int i)                                  { ltm.tm_min  += i; update(); }
			inline void        increment_hour(const int i)                                 { ltm.tm_hour += i; update(); }
			inline void        increment_day(const int i)                                  { ltm.tm_mday += i; update(); }
			inline void        increment_month(const int i)                                { ltm.tm_mon  += i; update(); }
			inline void        increment_year(const int i)                                 { ltm.tm_year += i; update(); }
			inline void        decrement_sec(const int i)                                  { ltm.tm_sec  -= i; update(); }
			inline void        decrement_min(const int i)                                  { ltm.tm_min  -= i; update(); }
			inline void        decrement_hour(const int i)                                 { ltm.tm_hour -= i; update(); }
			inline void        decrement_day(const int i)                                  { ltm.tm_mday -= i; update(); }
			inline void        decrement_month(const int i)                                { ltm.tm_mon  -= i; update(); }
			inline void        decrement_year(const int i)                                 { ltm.tm_year -= i; update(); }

			inline TimeObject &operator=(const tm &in)                                     { set_tm(in);     return (*this); }
			inline TimeObject &operator=(const TimeObject &TO)                             { set_tm(TO.ltm); return (*this); }
			inline time_t      operator+(const TimeObject &TO) const                       { return t + TO.t; }
			inline time_t      operator-(const TimeObject &TO) const                       { return t - TO.t; }

			inline bool        operator==(const TimeObject &TO) const                      { return t == TO.t; }
			inline bool        operator< (const TimeObject &TO) const                      { return t <  TO.t; }
			inline bool        operator<=(const TimeObject &TO) const                      { return t <= TO.t; }
			inline bool        operator> (const TimeObject &TO) const                      { return t >  TO.t; }
			inline bool        operator>=(const TimeObject &TO) const                      { return t >= TO.t; }

		private:
			time_t  t;
			std::tm ltm;

			inline void        now()                                                       { t = time(0);      ltm = *localtime(&t); ltm.tm_isdst = 0; }
			inline void        update()                                                    { t = timegm(&ltm); ltm = *gmtime(&t);    ltm.tm_isdst = 0; }
	};

	class TimeObjectList {

		public:
			inline TimeObjectList()                                {}
			inline TimeObjectList(const int nn)                    { n = nn; allocate(); }
			inline TimeObjectList(const int nn, TimeObject *TO)    { n = nn; allocate(); fill(TO); }
			inline TimeObjectList(const TimeObjectList &TOLi)      { n = TOLi.n; allocate(); fill(TOLi.TOList); }
			inline ~TimeObjectList()                               { clear(); }

			inline TimeObjectList(const TimeObject &start, const TimeObject &end, const char *delta);

			inline int            len() const                      { return n; }
			inline void           clear()                          { deallocate(); n = 0; }
			inline TimeObject    *data()                           { return TOList; }
			inline bool           isempty()                        { return n==0; }
			inline std::string    as_string(const char *fmt)       { std::string str = std::string("[")+TOList[0].as_string(fmt); for(int ii=1;ii<n;++ii) str+=std::string(", ")+TOList[ii].as_string(fmt); str+=std::string("]"); return str; }
			inline void           sort()                           { std::sort(begin(),end()); }
			inline int            find(const TimeObject &TO);
			inline TimeObjectList merge(TimeObjectList &TOLi, bool exact);

			inline TimeObject  operator[](int i) const             { return (i>=0) ? TOList[i] : TOList[n+i]; }
			inline TimeObject &operator[](int i)                   { return (i>=0) ? TOList[i] : TOList[n+i]; }

			inline TimeObjectList &operator=(const TimeObjectList &TOLi) { n = TOLi.n; allocate(); fill(TOLi.TOList); return (*this); }

			class iterator {
				public:
					typedef TimeObject              value_type;
					typedef std::ptrdiff_t          difference_type;
					typedef TimeObject*             pointer;
					typedef TimeObject&             reference;
					typedef std::input_iterator_tag iterator_category;

					inline iterator() : TOL(nullptr), i(0) {}
					inline iterator(TimeObjectList* TOLi, int ii) : TOL(TOLi), i(ii) {}

					inline       TimeObject& operator*()            { return (*TOL)[i]; }
					inline const TimeObject& operator*() const      { return (*TOL)[i]; }

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
					TimeObjectList* TOL;
					int  i;
			};

			inline iterator begin()                               { return iterator{this,0}; }
			inline iterator end()                                 { return iterator{this,n}; }

		private:
			int n;
			bool alloc;
			TimeObject *TOList;

			inline void allocate()           { if(~alloc) {TOList = new TimeObject[n]; alloc=true;} }
			inline void deallocate()         { if(alloc)  {delete [] TOList;           alloc=false;} }
			inline void fill(TimeObject *TO) { if(alloc) std::memcpy(TOList,TO,n*sizeof(TimeObject)); }
	};

	inline TimeObjectList::TimeObjectList(const TimeObject &start, const TimeObject &end, const char *delta) {
		/*
			Given a start and an end point, construct a time object list spaced every
			delta. Options for delta are:
				> "secs   = x" : every x seconds;
				> "mins   = x" : every x minutes;
				> "hours  = x" : every x hours;
				> "days   = x" : every x days;
				> "months = x" : every x months;
				> "years  = x" : every x years;
		*/
		// First parse the string
		std::string str(delta); size_t pos = str.find("=");
		std::string delta_type = str.substr(0,pos);
		int d = std::stoi(str.substr(pos+1,str.size()));

		TimeObject TO_next(start), TO_end(end);
		std::vector<TimeObject> TO_list;
		n = 0;

		do{
			TO_list.push_back( TimeObject(TO_next) ); n++;
			if (delta_type == std::string("secs"))
				TO_next.increment_sec(d);
			else if (delta_type == std::string("mins"))
				TO_next.increment_min(d);
			else if (delta_type == std::string("hours"))
				TO_next.increment_hour(d);
			else if (delta_type == std::string("days"))
				TO_next.increment_day(d);
			else if (delta_type == std::string("months"))
				TO_next.increment_month(d);
			else if (delta_type == std::string("years"))
				TO_next.increment_year(d);
		} while(TO_next <= TO_end);

		allocate();
		fill(&TO_list[0]);
	}

	inline int TimeObjectList::find(const TimeObject &TO) {
		/*
			Finds the nearest
			Argument:
				a datetime object
			Returns the index of the nearest element		
		*/
		int posmin = -1; time_t diffmin = 9999999;
		for (int ii=0;ii<n;++ii)  {
			time_t diff = TOList[ii] - TO;
			if (diff < diffmin) { diffmin = diff; posmin = ii; }
		}
		return posmin;
	}

	inline TimeObjectList TimeObjectList::merge(TimeObjectList &TOLi, bool exact=false) {
		/*
			Merge with a datetimelist
			Returns a datetimelist with all the dates
			from both lists without repetitions	
		*/
		std::vector<TimeObject> list1(n+TOLi.len()),list2(n+TOLi.len());
		std::memcpy(list1.data(),TOList,n*sizeof(TimeObject));                 // from 0 to n copy this->TOList
		std::memcpy(list1.data()+n,TOLi.data(),TOLi.len()*sizeof(TimeObject)); // from n to TOLi.len() copy TOLi
		std::sort(list1.begin(),list1.end());                                  // sort the list
		// Now remove duplicates from list1 to list2
		std::vector<TimeObject>::iterator iter;
		iter = (exact) ? std::unique_copy(list1.begin(),list1.end(),list2.begin()) :
			std::unique_copy(list1.begin(),list1.end(),list2.begin(), 
				[](TimeObject T1,TimeObject T2){ return abs(T1 - T2) < 24*3600; } );
		list2.resize( std::distance(list2.begin(),iter) );
		// List2 contains the unique values of the two lists
		return TimeObjectList((int)(list2.size()),list2.data());
	}
}

#endif