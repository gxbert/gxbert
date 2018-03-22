//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GXAllocator.hh 88444 2015-02-20 13:43:16Z gcosmo $
//
// 
// ------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// A class for fast allocation of objects to the heap through a pool of
// chunks organised as linked list. It's meant to be used by associating
// it to the object to be allocated and defining for it new and delete
// operators via MallocSingle() and FreeSingle() methods.
       
//      ---------------- GXAllocator ----------------
//
// Author: G.Cosmo (CERN), November 2000
// ------------------------------------------------------------

#ifndef GXAllocator_h
#define GXAllocator_h 1

#include <cstddef>
#include <typeinfo>

#include "GXAllocatorPool.hh"

class GXAllocatorBase
{
  public:
    GXAllocatorBase(); 
    virtual ~GXAllocatorBase();
    virtual void ResetStorage()=0;
    virtual size_t GetAllocatedSize() const=0;
    virtual int GetNoPages() const=0;
    virtual size_t GetPageSize() const=0;
    virtual void IncreasePageSize( unsigned int sz )=0;
    virtual const char* GetPoolType() const=0;
};

template <class Type>
class GXAllocator : public GXAllocatorBase
{
  public:  // with description

    GXAllocator() throw();
    ~GXAllocator() throw();
      // Constructor & destructor

    inline Type* MallocSingle();
    inline void FreeSingle(Type* anElement);
      // Malloc and Free methods to be used when overloading
      // new and delete operators in the client <Type> object

    inline void ResetStorage();
      // Returns allocated storage to the free store, resets allocator.
      // Note: contents in memory are lost using this call !

    inline size_t GetAllocatedSize() const;
      // Returns the size of the total memory allocated
    inline int GetNoPages() const;
      // Returns the total number of allocated pages
    inline size_t GetPageSize() const;
      // Returns the current size of a page
    inline void IncreasePageSize( unsigned int sz );
      // Resets allocator and increases default page size of a given factor

    inline const char* GetPoolType() const;
      // Returns the type_info Id of the allocated type in the pool

  public:  // without description

    // This public section includes standard methods and types
    // required if the allocator is to be used as alternative
    // allocator for STL containers.
    // NOTE: the code below is a trivial implementation to make
    //       this class an STL compliant allocator.
    //       It is anyhow NOT recommended to use this class as
    //       alternative allocator for STL containers !

    typedef Type value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef Type* pointer;
    typedef const Type* const_pointer;
    typedef Type& reference;
    typedef const Type& const_reference;

    template <class U> GXAllocator(const GXAllocator<U>& right) throw()
      : mem(right.mem) {}
      // Copy constructor

    pointer address(reference r) const { return &r; }
    const_pointer address(const_reference r) const { return &r; }
      // Returns the address of values

    pointer allocate(size_type n, void* = 0)
    {
      // Allocates space for n elements of type Type, but does not initialise
      //
      Type* mem_alloc = 0;
      if (n == 1)
        mem_alloc = MallocSingle();
      else
        mem_alloc = static_cast<Type*>(::operator new(n*sizeof(Type)));
      return mem_alloc;
    }
    void deallocate(pointer p, size_type n)
    {
      // Deallocates n elements of type Type, but doesn't destroy
      //
      if (n == 1)
        FreeSingle(p);
      else
        ::operator delete((void*)p);
      return;
    }

    void construct(pointer p, const Type& val) { new((void*)p) Type(val); }
      // Initialises *p by val
    void destroy(pointer p) { p->~Type(); }
      // Destroy *p but doesn't deallocate

    size_type max_size() const throw()
    {
      // Returns the maximum number of elements that can be allocated
      //
      return 2147483647/sizeof(Type);
    }

    template <class U>
    struct rebind { typedef GXAllocator<U> other; };
      // Rebind allocator to type U

    GXAllocatorPool mem;
      // Pool of elements of sizeof(Type)

  private:

    const char* tname;
      // Type name identifier
};

// ------------------------------------------------------------
// Inline implementation
// ------------------------------------------------------------

// Initialization of the static pool
//
// template <class Type> GXAllocatorPool GXAllocator<Type>::mem(sizeof(Type));

// ************************************************************
// GXAllocator constructor
// ************************************************************
//
template <class Type>
GXAllocator<Type>::GXAllocator() throw()
  : mem(sizeof(Type))
{
  tname = typeid(Type).name();
}

// ************************************************************
// GXAllocator destructor
// ************************************************************
//
template <class Type>
GXAllocator<Type>::~GXAllocator() throw()
{
}

// ************************************************************
// MallocSingle
// ************************************************************
//
template <class Type>
Type* GXAllocator<Type>::MallocSingle()
{
  return static_cast<Type*>(mem.Alloc());
}

// ************************************************************
// FreeSingle
// ************************************************************
//
template <class Type>
void GXAllocator<Type>::FreeSingle(Type* anElement)
{
  mem.Free(anElement);
  return;
}

// ************************************************************
// ResetStorage
// ************************************************************
//
template <class Type>
void GXAllocator<Type>::ResetStorage()
{
  // Clear all allocated storage and return it to the free store
  //
  mem.Reset();
  return;
}

// ************************************************************
// GetAllocatedSize
// ************************************************************
//
template <class Type>
size_t GXAllocator<Type>::GetAllocatedSize() const
{
  return mem.Size();
}

// ************************************************************
// GetNoPages
// ************************************************************
//
template <class Type>
int GXAllocator<Type>::GetNoPages() const
{
  return mem.GetNoPages();
}

// ************************************************************
// GetPageSize
// ************************************************************
//
template <class Type>
size_t GXAllocator<Type>::GetPageSize() const
{
  return mem.GetPageSize();
}

// ************************************************************
// IncreasePageSize
// ************************************************************
//
template <class Type>
void GXAllocator<Type>::IncreasePageSize( unsigned int sz )
{
  ResetStorage();
  mem.GrowPageSize(sz); 
}

// ************************************************************
// GetPoolType
// ************************************************************
//
template <class Type>
const char* GXAllocator<Type>::GetPoolType() const
{
  return tname;
}

// ************************************************************
// operator==
// ************************************************************
//
template <class T1, class T2>
bool operator== (const GXAllocator<T1>&, const GXAllocator<T2>&) throw()
{
  return true;
}

// ************************************************************
// operator!=
// ************************************************************
//
template <class T1, class T2>
bool operator!= (const GXAllocator<T1>&, const GXAllocator<T2>&) throw()
{
  return false;
}

#endif
