#include "QvisBinnerWindow.h"

#include <BinnerFilter.h>
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
// Method: QvisBinnerWindow::QvisBinnerWindow
//
// Purpose: 
//   Constructor
//
// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

QvisBinnerWindow::QvisBinnerWindow(const int type,
                         Binner *subj,
                         const char *caption,
                         const char *shortName,
                         QvisNotepadArea *notepad)
    : QvisOperatorWindow(type,subj, caption, shortName, notepad)
{
    atts = subj;
}


// ****************************************************************************
// Method: QvisBinnerWindow::~QvisBinnerWindow
//
// Purpose: 
//   Destructor
//
// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

QvisBinnerWindow::~QvisBinnerWindow()
{
}


// ****************************************************************************
// Method: QvisBinnerWindow::CreateWindowContents
//
// Purpose: 
//   Creates the widgets for the window.
//
// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisBinnerWindow::CreateWindowContents()
{
    QGridLayout *mainLayout = new QGridLayout(topLayout, 3,2,  10, "mainLayout");


    dim1Label = new QLabel("dim1", central, "dim1Label");
    mainLayout->addWidget(dim1Label,0,0);
    dim1 = new QLineEdit(central, "dim1");
    connect(dim1, SIGNAL(returnPressed()),
            this, SLOT(dim1ProcessText()));
    mainLayout->addWidget(dim1, 0,1);

    dim2Label = new QLabel("dim2", central, "dim2Label");
    mainLayout->addWidget(dim2Label,1,0);
    dim2 = new QLineEdit(central, "dim2");
    connect(dim2, SIGNAL(returnPressed()),
            this, SLOT(dim2ProcessText()));
    mainLayout->addWidget(dim2, 1,1);

    dim3Label = new QLabel("dim3", central, "dim3Label");
    mainLayout->addWidget(dim3Label,2,0);
    dim3 = new QLineEdit(central, "dim3");
    connect(dim3, SIGNAL(returnPressed()),
            this, SLOT(dim3ProcessText()));
    mainLayout->addWidget(dim3, 2,1);

}


// ****************************************************************************
// Method: QvisBinnerWindow::UpdateWindow
//
// Purpose: 
//   Updates the widgets in the window when the subject changes.
//
// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisBinnerWindow::UpdateWindow(bool doAll)
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
          case 0: //dim1
            temp.sprintf("%d", atts->GetDim1());
            dim1->setText(temp);
            break;
          case 1: //dim2
            temp.sprintf("%d", atts->GetDim2());
            dim2->setText(temp);
            break;
          case 2: //dim3
            temp.sprintf("%d", atts->GetDim3());
            dim3->setText(temp);
            break;
        }
    }
}


// ****************************************************************************
// Method: QvisBinnerWindow::GetCurrentValues
//
// Purpose: 
//   Gets values from certain widgets and stores them in the subject.
//
// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//   
// ****************************************************************************

void
QvisBinnerWindow::GetCurrentValues(int which_widget)
{
    bool okay, doAll = (which_widget == -1);
    QString msg, temp;

    // Do dim1
    if(which_widget == 0 || doAll)
    {
        temp = dim1->displayText().simplifyWhiteSpace();
        okay = !temp.isEmpty();
        if(okay)
        {
            int val = temp.toInt(&okay);
            atts->SetDim1(val);
        }

        if(!okay)
        {
            msg.sprintf("The value of dim1 was invalid. "
                "Resetting to the last good value of %d.",
                atts->GetDim1());
            Message(msg);
            atts->SetDim1(atts->GetDim1());
        }
    }

    // Do dim2
    if(which_widget == 1 || doAll)
    {
        temp = dim2->displayText().simplifyWhiteSpace();
        okay = !temp.isEmpty();
        if(okay)
        {
            int val = temp.toInt(&okay);
            atts->SetDim2(val);
        }

        if(!okay)
        {
            msg.sprintf("The value of dim2 was invalid. "
                "Resetting to the last good value of %d.",
                atts->GetDim2());
            Message(msg);
            atts->SetDim2(atts->GetDim2());
        }
    }

    // Do dim3
    if(which_widget == 2 || doAll)
    {
        temp = dim3->displayText().simplifyWhiteSpace();
        okay = !temp.isEmpty();
        if(okay)
        {
            int val = temp.toInt(&okay);
            atts->SetDim3(val);
        }

        if(!okay)
        {
            msg.sprintf("The value of dim3 was invalid. "
                "Resetting to the last good value of %d.",
                atts->GetDim3());
            Message(msg);
            atts->SetDim3(atts->GetDim3());
        }
    }

}


//
// Qt Slot functions
//


void
QvisBinnerWindow::dim1ProcessText()
{
    GetCurrentValues(0);
    Apply();
}


void
QvisBinnerWindow::dim2ProcessText()
{
    GetCurrentValues(1);
    Apply();
}


void
QvisBinnerWindow::dim3ProcessText()
{
    GetCurrentValues(2);
    Apply();
}


