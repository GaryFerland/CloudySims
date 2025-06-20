/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MODULE_H_
#define MODULE_H_

class module;
class t_warnings;

/**
 * Singleton container that tracks all `module` instances.
 * 
 * Responsibilities:
 *  • Stores registered modules in a vector (`m_l`).
 *  • Provides methods to initialize (`zero()`) and document (`comment()`) each module.
 *  • Automatically populated via the base class constructor `module::module()`.
 */
class module_list : public Singleton<module_list>
{
	friend class Singleton<module_list>;
public:
	vector<module *> m_l;
protected:
	module_list() {}
public:
	void use(module *m)
	{
		m_l.push_back(m);
	}
	void zero() const;
	void comment(t_warnings&) const;
};

/** abstract base class that declares three pure virtual functions:
	•	zero() — for resetting state.
	•	comment(t_warnings&) — for writing module-specific output.
	•	chName() — returns the module’s name as a string.
	•	Registers itself with a global module_list on construction:*/
class module
{
public:
	module()
	{
		module_list::Inst().use(this);
	}
	virtual void zero() = 0;
	virtual void comment(t_warnings&) = 0;
	virtual const char* chName() const = 0;
	virtual ~module() {}
};

#endif /* MODULE_ */
