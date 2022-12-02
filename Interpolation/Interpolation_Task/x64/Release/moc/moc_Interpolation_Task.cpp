/****************************************************************************
** Meta object code from reading C++ file 'Interpolation_Task.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../Interpolation_Task.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Interpolation_Task.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Interpolation_Task_t {
    QByteArrayData data[8];
    char stringdata0[159];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Interpolation_Task_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Interpolation_Task_t qt_meta_stringdata_Interpolation_Task = {
    {
QT_MOC_LITERAL(0, 0, 18), // "Interpolation_Task"
QT_MOC_LITERAL(1, 19, 8), // "makePlot"
QT_MOC_LITERAL(2, 28, 0), // ""
QT_MOC_LITERAL(3, 29, 28), // "on_actionLoad_File_triggered"
QT_MOC_LITERAL(4, 58, 39), // "on_actionGenerate_random_data..."
QT_MOC_LITERAL(5, 98, 31), // "on_comboBox_currentIndexChanged"
QT_MOC_LITERAL(6, 130, 4), // "arg1"
QT_MOC_LITERAL(7, 135, 23) // "on_actionInfo_triggered"

    },
    "Interpolation_Task\0makePlot\0\0"
    "on_actionLoad_File_triggered\0"
    "on_actionGenerate_random_data_triggered\0"
    "on_comboBox_currentIndexChanged\0arg1\0"
    "on_actionInfo_triggered"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Interpolation_Task[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x08 /* Private */,
       3,    0,   40,    2, 0x08 /* Private */,
       4,    0,   41,    2, 0x08 /* Private */,
       5,    1,   42,    2, 0x08 /* Private */,
       7,    0,   45,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    6,
    QMetaType::Void,

       0        // eod
};

void Interpolation_Task::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Interpolation_Task *_t = static_cast<Interpolation_Task *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->makePlot(); break;
        case 1: _t->on_actionLoad_File_triggered(); break;
        case 2: _t->on_actionGenerate_random_data_triggered(); break;
        case 3: _t->on_comboBox_currentIndexChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->on_actionInfo_triggered(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Interpolation_Task::staticMetaObject = { {
    &QMainWindow::staticMetaObject,
    qt_meta_stringdata_Interpolation_Task.data,
    qt_meta_data_Interpolation_Task,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Interpolation_Task::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Interpolation_Task::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Interpolation_Task.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int Interpolation_Task::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
