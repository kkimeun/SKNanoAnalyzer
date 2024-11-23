#ifndef Muon_h
#define Muon_h

#include "TString.h"
#include "Lepton.h"

// Need update
// - TuneP object
// - momentum scale
// - rochester correction
// - chi2
// - truth matching info
// - Analysis dependent IDs

class Muon: public Lepton {
public:
    Muon();
    ~Muon();

    // Boolean IDs
    enum class BooleanID {NONE, LOOSE, MEDIUM, MEDIUMPROMPT, TIGHT, SOFT, SOFTMVA, TRIGGERLOOSE};

    void SetBIDBit(BooleanID id, bool idbit);
    inline bool isPOGTightId() const {return j_tightId;}
    inline bool isPOGMediumId() const {return j_mediumId;}
    inline bool isPOGMediumPromptId() const {return j_mediumPromptId;}
    inline bool isPOGLooseId() const {return j_looseId;}
    inline bool isPOGSoftId() const {return j_softId;}
    inline bool isPOGSoftMvaId() const {return j_softMvaId;}
    inline bool isPOGTriggerIdLoose() const {return j_triggerIdLoose;}

    // Unsigned char IDs
    enum class WorkingPointID {NONE, HIGHPT, MINIISO, MULTIISO, MVAMU, PFISO, PUPPIISO, TKISO};
    enum class WorkingPoint {NONE, VLOOSE, LOOSE, MEDIUM, TIGHT, VTIGHT, VVTIGHT};

    enum class MuonID
    {
        NOCUT,
        POG_TIGHT,
        POG_MEDIUM,
        POG_MEDIUM_PROMPT,
        POG_LOOSE,
        POG_SOFT,
        POG_SOFT_MVA,
        POG_TRIGGER_LOOSE,
        POG_TRACKER_HIGH_PT,
        POG_GLOBAL_HIGH_PT,
        POG_MINISO_LOOSE,
        POG_MINISO_MEDIUM,
        POG_MINISO_TIGHT,
        POG_MINISO_VTIGHT,
        POG_MULTISO_LOOSE,
        POG_MULTISO_MEDIUM,
        POG_MVA_MU_MEDIUM,
        POG_MVA_MU_TIGHT,
        POG_PFISO_VLOOSE,
        POG_PFISO_LOOSE,
        POG_PFISO_MEDIUM,
        POG_PFISO_TIGHT,
        POG_PFISO_VTIGHT,
        POG_PFISO_VVTIGHT,
        POG_PUPPIISO_LOOSE,
        POG_PUPPIISO_MEDIUM,
        POG_PUPPIISO_TIGHT,
        POG_TKISO_LOOSE,
        POG_TKISO_TIGHT
    };

    void SetWIDBit(WorkingPointID id, unsigned char value);
    inline WorkingPoint HighPtId() const {return (WorkingPoint)j_highPtId;}
    inline WorkingPoint MiniIsoId() const {return (WorkingPoint)j_miniIsoId;}
    inline WorkingPoint MultiIsoId() const {return (WorkingPoint)j_multiIsoId;}
    inline WorkingPoint MvaMuId() const {return (WorkingPoint)j_mvaMuId;}
    //inline WorkingPoint MvaLowPtId() const {return (WorkingPoint)j_mvaLowPtId;}
    inline WorkingPoint PfIsoId() const {return (WorkingPoint)j_pfIsoId;}
    inline WorkingPoint PuppiIsoId() const {return (WorkingPoint)j_puppiIsoId;}
    inline WorkingPoint TkIsoId() const {return (WorkingPoint)j_tkIsoId;}

    // MVA ID scores
    enum class MVAID {NONE, SOFTMVA, MVALOWPT, MVATTH};

    void SetMVAID(MVAID id, float score);
    inline float SoftMva() const {return j_softMva;}
    inline float MvaLowPt() const {return j_mvaLowPt;}
    inline float MvaTTH() const {return j_mvaTTH;}

    // ID helper functions
    bool PassID(const MuonID ID) const;
    bool PassID(const TString ID) const;

private:
    bool j_isTracker, j_isStandalone, j_isGlobal;
    bool j_looseId, j_mediumId, j_mediumPromptId, j_tightId, j_softId, j_softMvaId, j_triggerIdLoose;
    unsigned char j_highPtId, j_miniIsoId, j_multiIsoId, j_mvaMuId, j_pfIsoId, j_puppiIsoId, j_tkIsoId;
    float j_softMva, j_mvaLowPt, j_mvaTTH;

    ClassDef(Muon, 1);
};

#endif
