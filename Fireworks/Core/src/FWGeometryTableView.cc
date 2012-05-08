// -*- C++ -*-
//
// Package:     Core
// Class  :     FWGeometryTableView
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  
//         Created:  Wed Jan  4 00:05:34 CET 2012
// $Id: FWGeometryTableView.cc,v 1.34 2012/05/04 03:00:38 amraktad Exp $
//

// system include files
#include <boost/bind.hpp>

// user include files
#include "Fireworks/Core/src/FWGeometryTableView.h"
#include "Fireworks/Core/src/FWGeoTopNodeScene.h"
#include "Fireworks/Core/interface/FWGeometryTableViewManager.h"
#include "Fireworks/Core/interface/FWViewType.h"
#include "Fireworks/Core/interface/FWGeometryTableManagerBase.h"
#include "Fireworks/Core/interface/CmsShowViewPopup.h"
#include "Fireworks/Core/src/FWGeometryTableManager.h"
#include "Fireworks/Core/interface/fwLog.h"

#include "Fireworks/Core/src/FWGUIValidatingTextEntry.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/src/FWValidatorBase.h"
#include "Fireworks/Core/src/FWEveDetectorGeo.h"

#include "KeySymbols.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGListBox.h"
#include "TGLViewer.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"

#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEveSceneInfo.h"
#include "TEveManager.h"
#include "TGeoManager.h"
#include "TGLCamera.h"


//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================


class FWGeoMaterialValidator : public FWValidatorBase 
{
public:
   struct Material
   {
      TGeoMaterial* g;
      std::string n;
      bool operator< (const Material& x) const { return n < x.n ;}
      Material( TGeoMaterial* x) {  g= x; n = x ? x->GetName() : "<show-all>";}
   };

   FWGeometryTableView* m_browser;
   mutable std::vector<Material> m_list;

   FWGeoMaterialValidator( FWGeometryTableView* v) { m_browser = v;}
   virtual ~FWGeoMaterialValidator() {}


   virtual void fillOptions(const char* iBegin, const char* iEnd, std::vector<std::pair<boost::shared_ptr<std::string>, std::string> >& oOptions) const 
   {
      oOptions.clear();
      std::string part(iBegin,iEnd);
      unsigned int part_size = part.size();

      m_list.clear();
      m_list.push_back(0);

      FWGeometryTableManagerBase::Entries_i it = m_browser->getTableManager()->refEntries().begin();
      std::advance(it, m_browser->getTopNodeIdx());
      int nLevel = it->m_level;
      it++;
      while (it->m_level > nLevel)
      {
         TGeoMaterial* g = it->m_node->GetVolume()->GetMaterial();
         bool duplicate = false;
         for (std::vector<Material>::iterator j = m_list.begin(); j!=m_list.end(); ++j) {
            if (j->g == g) {
               duplicate = true;
               break;
            }
         }
         if (!duplicate)
            m_list.push_back(g);

         ++it;
      }
      std::vector<Material>::iterator startIt = m_list.begin();
      startIt++;
      std::sort(startIt, m_list.end());

      std::string h = "";
      oOptions.push_back(std::make_pair(boost::shared_ptr<std::string>(new std::string(m_list.begin()->n)), h));
      for (std::vector<Material>::iterator i = startIt; i!=m_list.end(); ++i)
      {
         if (part == (*i).n.substr(0,part_size))
         {
            //  std::cout << i->n <<std::endl;
            oOptions.push_back(std::make_pair(boost::shared_ptr<std::string>(new std::string((*i).n)), (*i).n.substr(part_size, (*i).n.size()-part_size)));
         }
      }
   }

   bool isStringValid(std::string& exp) 
   {
      if (exp.empty()) return true;

      for (std::vector<Material>::iterator i = m_list.begin(); i != m_list.end(); ++i)
      {
         if (exp == (*i).n) 
            return true;
      }
      return false;
   }
};


//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================

FWGeometryTableView::FWGeometryTableView(TEveWindowSlot* iParent, FWColorManager* colMng)
   : FWGeometryTableViewBase(iParent, FWViewType::kGeometryTable, colMng),
     m_tableManager(0),
     m_filterEntry(0),
     m_filterValidator(0),
     m_mode(this, "Mode", 0l, 0l, 1l),
     m_disableTopNode(this,"HideTopNode", true),
     m_visLevel(this,"VisLevel", 3l, 1l, 100l),
     m_filter(this,"Materials", std::string()),
     m_filterType(this,"FilterType", 0l, 0l, 3l),
     m_visLevelFilter(this,"IgnoreVisLevelOnFilter", true),
     m_selectRegion(this, "SelectNearCameraCenter", false),
     m_regionRadius(this, "SphereRadius", 10.0, 1.0, 300.0),
     m_proximityAlgo(this, "Proximity algorithm", 1l, 0l, 1l)
{
   FWGeoTopNodeGLScene *gls = new FWGeoTopNodeGLScene(0);
#if ROOT_VERSION_CODE < ROOT_VERSION(5,32,0)
   m_eveScene  = new  FWGeoTopNodeEveScene(gls, "TopGeoNodeScene", "");
#else
   m_eveScene  = new  TEveScene(gls, "TopGeoNodeScene", "");
#endif
   gEve->GetScenes()->AddElement(m_eveScene);

   m_eveTopNode = new  FWEveDetectorGeo(this);
   m_eveTopNode->IncDenyDestroy();
   m_eveTopNode->SetPickable(true);
   m_eveScene->AddElement(m_eveTopNode);

   gls->m_eveTopNode = m_eveTopNode;
   m_eveTopNode->m_scene   = gls;

   // top row
   TGHorizontalFrame *hp = new TGHorizontalFrame(m_frame);
   {
      TGTextButton *rb = new TGTextButton (hp, "CdTop");
      hp->AddFrame(rb, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0) );
      rb->Connect("Clicked()","FWGeometryTableViewBase",this,"cdTop()");
   } 
   {
      TGTextButton *rb = new TGTextButton (hp, "CdUp");
      hp->AddFrame(rb, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
      rb->Connect("Clicked()","FWGeometryTableViewBase",this,"cdUp()");
   }
   {
      m_viewBox = new FWViewCombo(hp, this);
      hp->AddFrame( m_viewBox,new TGLayoutHints(kLHintsExpandY, 2, 2, 0, 0));
   }
   {
      hp->AddFrame(new TGLabel(hp, "Filter:"), new TGLayoutHints(kLHintsBottom, 10, 0, 0, 2));
      m_filterEntry = new FWGUIValidatingTextEntry(hp);
      m_filterEntry->SetHeight(20);
      m_filterValidator = new FWGeoMaterialValidator(this);
      m_filterEntry->setValidator(m_filterValidator);
      hp->AddFrame(m_filterEntry, new TGLayoutHints(kLHintsExpandX,  1, 2, 1, 0));
      m_filterEntry->setMaxListBoxHeight(150);
      m_filterEntry->getListBox()->Connect("Selected(int)", "FWGeometryTableView",  this, "filterListCallback()");
      m_filterEntry->Connect("ReturnPressed()", "FWGeometryTableView",  this, "filterTextEntryCallback()");

      gVirtualX->GrabKey( m_filterEntry->GetId(),gVirtualX->KeysymToKeycode((int)kKey_A),  kKeyControlMask, true);
   }
   m_frame->AddFrame(hp,new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 4, 2, 2, 0));

   m_tableManager = new FWGeometryTableManager(this);
   {
      TEveGeoManagerHolder gmgr( FWGeometryTableViewManager::getGeoMangeur());
      m_tableManager->loadGeometry( gGeoManager->GetTopNode(), gGeoManager->GetListOfVolumes());
   }
   cdTop();

   m_mode.addEntry(kNode,   "Node");
   m_mode.addEntry(kVolume, "Volume");

   m_mode.changed_.connect(boost::bind(&FWGeometryTableView::refreshTable3D,this));
   m_autoExpand.changed_.connect(boost::bind(&FWGeometryTableView::autoExpandCallback, this));
   m_visLevel.changed_.connect(boost::bind(&FWGeometryTableView::refreshTable3D,this));
   
   
   
   m_filterType.addEntry(kFilterMaterialName,   "MaterialName");
   m_filterType.addEntry(kFilterMaterialTitle,  "MaterialTitle");
   m_filterType.addEntry(kFilterShapeName,      "ShapeName");
   m_filterType.addEntry(kFilterShapeClassName, "ShapeClassName");

   m_visLevelFilter.changed_.connect(boost::bind(&FWGeometryTableView::refreshTable3D,this));

   m_disableTopNode.changed_.connect(boost::bind(&FWGeometryTableView::updateVisibilityTopNode,this));
   postConst();

   m_proximityAlgo.addEntry(kBBoxCenter,  "BBox center");
   m_proximityAlgo.addEntry(kBBoxSurface, "BBox surface");

   m_selectRegion.changed_.connect(boost::bind(&FWGeometryTableView::checkRegionOfInterest,this));
   m_regionRadius.changed_.connect(boost::bind(&FWGeometryTableView::checkRegionOfInterest,this));
   m_proximityAlgo.changed_.connect(boost::bind(&FWGeometryTableView::checkRegionOfInterest,this));
   
}

FWGeometryTableView::~FWGeometryTableView()
{}

//______________________________________________________________________________
void FWGeometryTableView::setPath(int parentIdx, std::string&)
{
   m_eveTopNode->clearSelection();

   m_topNodeIdx.set(parentIdx);
   getTableManager()->refEntries().at(getTopNodeIdx()).setBitVal(FWGeometryTableManagerBase::kVisNodeSelf,!m_disableTopNode.value() );
   getTableManager()->setLevelOffset(getTableManager()->refEntries().at(getTopNodeIdx()).m_level);

   checkExpandLevel();
   refreshTable3D(); 
}

//______________________________________________________________________________
FWGeometryTableManagerBase* FWGeometryTableView::getTableManager()
{
   return m_tableManager;
}

//______________________________________________________________________________
void FWGeometryTableView::autoExpandCallback()
{ 
   if (!m_enableRedraw) return;
   checkExpandLevel();
   getTableManager()->redrawTable(true);
}

//______________________________________________________________________________
void FWGeometryTableView::filterTextEntryCallback()
{
   // std::cout << "text entry click ed \n" ;
   std::string exp = m_filterEntry->GetText();
   updateFilter(exp);
}

//______________________________________________________________________________
void FWGeometryTableView::filterListCallback()
{ 
   // std::cout << "list click ed [" << m_filterEntry->GetText() << "] \n" ;

   std::string exp = m_filterEntry->GetText();
   updateFilter(exp);
}

//______________________________________________________________________________
void FWGeometryTableView::updateFilter(std::string& exp)
{
   // std::cout << "=FWGeometryTableViewBase::updateFilter()" << m_filterEntry->GetText() <<std::endl;
  
   
   if (exp.empty())
   {
      // std::cout << "FITLER OFF \n";
      for (FWGeometryTableManagerBase::Entries_i i = m_tableManager->refEntries().begin(); i !=  m_tableManager->refEntries().end(); ++i)
      {
         m_tableManager->setVisibility(*i, true);
         m_tableManager->setVisibilityChld(*i, true);
      }

      // NOTE: entry should be cleared automatically
      m_filterEntry->Clear();
   }
  
   m_filter.set(exp);
   m_tableManager->updateFilter(m_filterType.value());
   refreshTable3D();

}

//==============================================================================

void FWGeometryTableView::populateController(ViewerParameterGUI& gui) const
{
   gui.requestTab("Style").
      addParam(&m_disableTopNode).
      addParam(&m_mode).
      addParam(&m_autoExpand).
      addParam(&m_visLevel).
      separator().   
      addParam(&m_filterType).
      addParam(&m_visLevelFilter).
      separator().   
      addParam(&m_selectRegion).
      addParam(&m_regionRadius).
      addParam(&m_proximityAlgo);

      // addParam(&m_enableHighlight);
   
   
   FWGeometryTableViewBase::populateController(gui);
}


//------------------------------------------------------------------------------

/*
void FWGeometryTableView::setPath(int parentIdx, std::string& path)
{
   //   printf("Set Path to [%s], current node \n", path.c_str());
   m_topNodeIdx.set(parentIdx);
   getTableManager()->refEntries().at(getTopNodeIdx()).setBitVal(FWGeometryTableManagerBase::kVisNodeSelf,!m_disableTopNode.value() );
   getTableManager()->setLevelOffset(getTableManager()->refEntries().at(getTopNodeIdx()).m_level);

   m_eveTopNode->clearSelection(); 

   checkExpandLevel();
   refreshTable3D(); 
   FWGUIManager::getGUIManager()->updateStatus(path.c_str());
   }*/

//--------------------------------------------------------------
bool viewIsChecked(TEveViewer* v, TEveElement* el)
{
   if (strstr( v->GetElementName(), "3D") )
   {
      for (TEveElement::List_i eit = v->BeginChildren(); eit != v->EndChildren(); ++eit )
      {
         TEveScene* s = ((TEveSceneInfo*)*eit)->GetScene();
         if (el && s->HasChildren() && s->FirstChild() == el) 
            return true;
      }
      
   }
   return false;
}

void FWGeometryTableView::checkRegionOfInterest()
{
   if (m_selectRegion.value())
   {
      double* center = 0;
      for (TEveElement::List_i it = gEve->GetViewers()->BeginChildren(); it != gEve->GetViewers()->EndChildren(); ++it)
      { 
         TEveViewer* v = ((TEveViewer*)(*it));
         if (viewIsChecked(v, m_eveTopNode))
         {
            if (center) {
               fwLog(fwlog::kWarning) << "Center picked from first view \n";
            } else {
               center = v->GetGLViewer()->CurrentCamera().GetCenterVec();
               fwLog(fwlog::kInfo) << Form("Center picked (%.1f, %.1f, %.1f) from first selected 3D view \n", 
                                           center[0], center[1], center[2]);
            }
         }
      } 

      if (! center)
      {
         fwLog(fwlog::kError) << "No 3D view selected \n";
         return;
      }
      
      m_tableManager->checkRegionOfInterest(center, m_regionRadius.value(), m_proximityAlgo.value());
   }
   else 
   {
      m_tableManager->resetRegionOfInterest();
   }

   refreshTable3D();
}
//------------------------------------------------------------------------------

void FWGeometryTableView::setFrom(const FWConfiguration& iFrom)
{ 
   m_enableRedraw = false;
   for (const_iterator it =begin(), itEnd = end(); it != itEnd; ++it)
   {
      //      printf("set from %s \n",(*it)->name().c_str() );
      (*it)->setFrom(iFrom);
   }  
   m_viewersConfig = iFrom.valueForKey("Viewers");

   cdNode(m_topNodeIdx.value());
   m_enableRedraw = true;

   checkExpandLevel();
   refreshTable3D();
   /*
   getTableManager()->redrawTable();
   m_eveTopNode->ElementChanged();
   gEve->FullRedraw3D(false, true);
   */
}

//------------------------------------------------------------------------------

void FWGeometryTableView::updateVisibilityTopNode()
{
   getTableManager()->refEntries().at(getTopNodeIdx()).setBitVal(FWGeometryTableManagerBase::kVisNodeSelf,!m_disableTopNode.value() );
   refreshTable3D();
}
