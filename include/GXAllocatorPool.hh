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
// $Id: GXAllocatorPool.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// 
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
// Class description:
//
// Class implementing a memory pool for fast allocation and deallocation
// of memory chunks.  The size of the chunks for small allocated objects
// is fixed to 1Kb and takes into account of memory alignment; for large
// objects it is set to 10 times the object's size.
// The implementation is derived from: B.Stroustrup, The C++ Programming
// Language, Third Edition.

//           -------------- GXAllocatorPool ----------------
//
// Author: G.Cosmo (CERN), November 2000
// -------------------------------------------------------------------

#ifndef GXAllocatorPool_h
#define GXAllocatorPool_h 1

class GXAllocatorPool
{
  public:

    explicit GXAllocatorPool( unsigned int n=0 );
      // Create a pool of elements of size n
    ~GXAllocatorPool();
      // Destructor. Return storage to the free store

    inline void* Alloc();
      // Allocate one element
    inline void  Free( void* b );
      // Return an element back to the pool

    inline unsigned int  Size() const;
      // Return storage size
    void  Reset();
      // Return storage to the free store

    inline int  GetNoPages() const;
      // Return the total number of allocated pages
    inline unsigned int  GetPageSize() const;
      // Accessor for default page size
    inline void GrowPageSize( unsigned int factor );
      // Increase default page size by a given factor

  private:

    GXAllocatorPool(const GXAllocatorPool& right);
      // Provate copy constructor
    GXAllocatorPool& operator= (const GXAllocatorPool& right);
      // Private equality operator

    struct G4PoolLink
    {
      G4PoolLink* next;
    };
    class G4PoolChunk
    {
      public:
        explicit G4PoolChunk(unsigned int sz)
          : size(sz), mem(new char[size]), next(0) {;}
        ~G4PoolChunk() { delete [] mem; }
        const unsigned int size;
        char* mem;
        G4PoolChunk* next;
    };

    void Grow();
      // Make pool larger

  private:

    const unsigned int esize;
    unsigned int csize;
    G4PoolChunk* chunks;
    G4PoolLink* head;
    int nchunks;
};

// ------------------------------------------------------------
// Inline implementation
// ------------------------------------------------------------

// ************************************************************
// Alloc
// ************************************************************
//
inline void*
GXAllocatorPool::Alloc()
{
  if (head==0) { Grow(); }
  G4PoolLink* p = head;  // return first element
  head = p->next;
  return p;
}

// ************************************************************
// Free
// ************************************************************
//
inline void
GXAllocatorPool::Free( void* b )
{
  G4PoolLink* p = static_cast<G4PoolLink*>(b);
  p->next = head;        // put b back as first element
  head = p;
}

// ************************************************************
// Size
// ************************************************************
//
inline unsigned int
GXAllocatorPool::Size() const
{
  return nchunks*csize;
}

// ************************************************************
// GetNoPages
// ************************************************************
//
inline int
GXAllocatorPool::GetNoPages() const
{
  return nchunks;
}

// ************************************************************
// GetPageSize
// ************************************************************
//
inline unsigned int
GXAllocatorPool::GetPageSize() const
{
  return csize;
}

// ************************************************************
// GrowPageSize
// ************************************************************
//
inline void
GXAllocatorPool::GrowPageSize( unsigned int sz )
{
  csize = (sz) ? sz*csize : csize; 
}

#endif
