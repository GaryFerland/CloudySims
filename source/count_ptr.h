/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#ifndef COUNT_PTR_H_
#define COUNT_PTR_H_

#include <algorithm>

// Temporary stand-in for shared_ptr<> until we can rely on it being
// available as standard
template <class T>
class count_ptr {
private:
	T* m_ptr;
	long *m_count;
public:
	// Constructors
	explicit count_ptr(T* ptr = 0)
		: m_ptr(ptr), m_count(new long(1))
		{	}
	count_ptr(const count_ptr<T> &p)
		: m_ptr(p.m_ptr), m_count(p.m_count)
		{
			++*m_count;
		}
	~count_ptr() throw ()
		{
			cancel();
		}
	// Modifiers
	count_ptr<T>& operator=(count_ptr<T> temp)
		{
			swap(temp);
			return *this;
		}
	// Don't implement arithmetic -- would corrupt the owned pointer
	// Accessors
	T& operator*() const
		{
			return *m_ptr;
		}
	T* operator->() const throw()
		{
			return m_ptr;
		}
	int equiv(const count_ptr<T>& p) const throw()
		{
			return m_ptr == p.m_ptr;
		}
	T* get_ptr() const throw()
		{
			return m_ptr;
		}
	long count() const
		{
			return *m_count;
		}
	void swap (count_ptr<T> &a) throw()
		{
			std::swap(m_ptr,a.m_ptr);
			std::swap(m_count,a.m_count);
		}
	int compare (const count_ptr<T> &a) const
		{
			if (m_ptr < a.m_ptr)
				return -1;
			else if (m_ptr > a.m_ptr)
				return 1;
			return 0;
		}
private:
	void cancel () 
		{
			if (0 == --*m_count)
			{
				delete m_count;
				delete m_ptr;
			}
		}
};

template <class T>
inline void swap (count_ptr<T> &a, count_ptr<T> &b)
{
	a.swap(b);
}

template <class T>
inline bool operator< (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return a.compare(b) < 0;
}
template <class T>
inline bool operator> (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return a.compare(b) > 0;
}
template <class T>
inline bool operator<= (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return a.compare(b) <= 0;
}
template <class T>
inline bool operator>= (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return a.compare(b) >= 0;
}
template <class T>
inline bool operator== (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return a.compare(b) == 0;
}
template <class T>
inline bool operator!= (const count_ptr<T> &a, const count_ptr<T> &b)
{
	return !(a == b);
}

#endif /* COUNT_PTR_H_ */
