// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooBernsteinSidebandDictionary

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "RooBernsteinSideband.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RooBernsteinSideband(void *p = 0);
   static void *newArray_RooBernsteinSideband(Long_t size, void *p);
   static void delete_RooBernsteinSideband(void *p);
   static void deleteArray_RooBernsteinSideband(void *p);
   static void destruct_RooBernsteinSideband(void *p);
   static void streamer_RooBernsteinSideband(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooBernsteinSideband*)
   {
      ::RooBernsteinSideband *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooBernsteinSideband >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooBernsteinSideband", ::RooBernsteinSideband::Class_Version(), "RooBernsteinSideband.h", 20,
                  typeid(::RooBernsteinSideband), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooBernsteinSideband::Dictionary, isa_proxy, 16,
                  sizeof(::RooBernsteinSideband) );
      instance.SetNew(&new_RooBernsteinSideband);
      instance.SetNewArray(&newArray_RooBernsteinSideband);
      instance.SetDelete(&delete_RooBernsteinSideband);
      instance.SetDeleteArray(&deleteArray_RooBernsteinSideband);
      instance.SetDestructor(&destruct_RooBernsteinSideband);
      instance.SetStreamerFunc(&streamer_RooBernsteinSideband);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooBernsteinSideband*)
   {
      return GenerateInitInstanceLocal((::RooBernsteinSideband*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooBernsteinSideband*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooBernsteinSideband::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooBernsteinSideband::Class_Name()
{
   return "RooBernsteinSideband";
}

//______________________________________________________________________________
const char *RooBernsteinSideband::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinSideband*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooBernsteinSideband::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinSideband*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooBernsteinSideband::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinSideband*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooBernsteinSideband::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinSideband*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooBernsteinSideband::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooBernsteinSideband.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> _maxDegree1;
      R__b >> _maxDegree2;
      R__b >> _maxDegree3;
      R__b >> _maxDegreeV;
      R__b.CheckByteCount(R__s, R__c, RooBernsteinSideband::IsA());
   } else {
      R__c = R__b.WriteVersion(RooBernsteinSideband::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << _maxDegree1;
      R__b << _maxDegree2;
      R__b << _maxDegree3;
      R__b << _maxDegreeV;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooBernsteinSideband(void *p) {
      return  p ? new(p) ::RooBernsteinSideband : new ::RooBernsteinSideband;
   }
   static void *newArray_RooBernsteinSideband(Long_t nElements, void *p) {
      return p ? new(p) ::RooBernsteinSideband[nElements] : new ::RooBernsteinSideband[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooBernsteinSideband(void *p) {
      delete ((::RooBernsteinSideband*)p);
   }
   static void deleteArray_RooBernsteinSideband(void *p) {
      delete [] ((::RooBernsteinSideband*)p);
   }
   static void destruct_RooBernsteinSideband(void *p) {
      typedef ::RooBernsteinSideband current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooBernsteinSideband(TBuffer &buf, void *obj) {
      ((::RooBernsteinSideband*)obj)->::RooBernsteinSideband::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooBernsteinSideband

namespace {
  void TriggerDictionaryInitialization_RooBernsteinSidebandDictionary_Impl() {
    static const char* headers[] = {
"RooBernsteinSideband.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/root-6.08.06/include",
"/gwpool/users/dini/RooSideBandGit/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RooBernsteinSidebandDictionary dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooBernsteinSideband.h")))  RooBernsteinSideband;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RooBernsteinSidebandDictionary dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "RooBernsteinSideband.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooBernsteinSideband", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooBernsteinSidebandDictionary",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooBernsteinSidebandDictionary_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooBernsteinSidebandDictionary_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooBernsteinSidebandDictionary() {
  TriggerDictionaryInitialization_RooBernsteinSidebandDictionary_Impl();
}
