#include "QvisH5PartRandomSampleWindow.h"

#include <H5PartRandomSampleAttributes.h>
#include <ViewerProxy.h>

#include <qcheckbox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qspinbox.h>
#include <qvbox.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <QvisColorTableButton.h>
#include <QvisOpacitySlider.h>
#include <QvisColorButton.h>
#include <QvisLineStyleWidget.h>
#include <QvisLineWidthWidget.h>
#include <QvisVariableButton.h>

#include <stdio.h>
#include <string>

using std::string;

// ****************************************************************************
// Method: QvisH5PartRandomSampleWindow::QvisH5PartRandomSampleWindow
//
// Purpose: 
//   Constructor
//
// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

QvisH5PartRandomSampleWindow::QvisH5PartRandomSampleWindow(const int type,
                         H5PartRandomSampleAttributes *subj,
                         const char *caption,
                         const char *shortName,
                         QvisNotepadArea *notepad)
    : QvisOperatorWindow(type,subj, caption, shortName, notepad)
{
    atts = subj;
}


// ****************************************************************************
// Method: QvisH5PartRandomSampleWindow::~QvisH5PartRandomSampleWindow
//
// Purpose: 
//   Destructor
//
// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

QvisH5PartRandomSampleWindow::~QvisH5PartRandomSampleWindow()
{
}


// ****************************************************************************
// Method: QvisH5PartRandomSampleWindow::CreateWindowContents
//
// Purpose: 
//   Creates the widgets for the window.
//
// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisH5PartRandomSampleWindow::CreateWindowContents()
{
    QGridLayout *mainLayout = new QGridLayout(topLayout, 1,2,  10, "mainLayout");


    factorLabel = new QLabel("factor", central, "factorLabel");
    mainLayout->addWidget(factorLabel,0,0);
    factor = new QLineEdit(central, "factor");
    connect(factor, SIGNAL(returnPressed()),
            this, SLOT(factorProcessText()));
    mainLayout->addWidget(factor, 0,1);

}


// ****************************************************************************
// Method: QvisH5PartRandomSampleWindow::UpdateWindow
//
// Purpose: 
//   Updates the widgets in the window when the subject changes.
//
// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisH5PartRandomSampleWindow::UpdateWindow(bool doAll)
{
    QString temp;
    double r;

    for(int i = 0; i < atts->NumAttributes(); ++i)
    {
        if(!doAll)
        {
            if(!atts->IsSelected(i))
            {
                continue;
            }
        }

        const double         *dptr;
        const float          *fptr;
        const int            *iptr;
        const char           *cptr;
        const unsigned char  *uptr;
        const string         *sptr;
        QColor                tempcolor;
        switch(i)
        {
          case 0: //factor
            temp.setNum(atts->GetFactor());
            factor->setText(temp);
            break;
        }
    }
}


// ****************************************************************************
// Method: QvisH5PartRandomSampleWindow::GetCurrentValues
//
// Purpose: 
//   Gets values from certain widgets and stores them in the subject.
//
// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisH5PartRandomSampleWindow::GetCurrentValues(int which_widget)
{
    bool okay, doAll = (which_widget == -1);
    QString msg, temp;

    // Do factor
    if(which_widget == 0 || doAll)
    {
        temp = factor->displayText().simplifyWhiteSpace();
        okay = !temp.isEmpty();
        if(okay)
        {
            float val = temp.toFloat(&okay);
            atts->SetFactor(val);
        }

        if(!okay)
        {
            msg.sprintf("The value of factor was invalid. "
                "Resetting to the last good value of %g.",
                atts->GetFactor());
            Message(msg);
            atts->SetFactor(atts->GetFactor());
        }
    }

}


//
// Qt Slot functions
//


void
QvisH5PartRandomSampleWindow::factorProcessText()
{
    GetCurrentValues(0);
    Apply();
}


