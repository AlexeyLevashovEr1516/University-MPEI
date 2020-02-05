;******************************************************************************
;   This file is a basic template for creating relocatable assembly code for  *
;   a PIC18F2520. Copy this file into your project directory and modify or    *
;   add to it as needed.                                                      *
;                                                                             *
;   The PIC18FXXXX architecture allows two interrupt configurations. This     *
;   template code is written for priority interrupt levels and the IPEN bit   *
;   in the RCON register must be set to enable priority levels. If IPEN is    *
;   left in its default zero state, only the interrupt vector at 0x008 will   *
;   be used and the WREG_TEMP, BSR_TEMP and STATUS_TEMP variables will not    *
;   be needed.                                                                *
;                                                                             *
;   Refer to the MPASM User's Guide for additional information on the         *
;   features of the assembler and linker.                                     *
;                                                                             *
;   Refer to the PIC18Fxx20 Data Sheet for additional                         *
;   information on the architecture and instruction set.                      *
;                                                                             *
;******************************************************************************
;                                                                             *
;    Filename:                                                                *
;    Date:                                                                    *
;    File Version:                                                            *
;                                                                             *
;    Author:                                                                  *
;    Company:                                                                 *
;                                                                             * 
;******************************************************************************
;                                                                             *
;    Files required: P18F2520.INC                                             *
;                                                                             *
;******************************************************************************

	LIST P=18F2520, F=INHX32 ;directive to define processor and file format
	#include <p18f2520.INC>	 ;processor specific variable definitions

;******************************************************************************
;Configuration bits
;Microchip has changed the format for defining the configuration bits, please 
;see the .inc file for futher details on notation.  Below are a few examples.



;   Oscillator Selection:
    CONFIG	OSC = INTIO67, WDT = OFF, PBADEN = OFF;, MCLRE = ON, DEBUG = ON             ;LP

;******************************************************************************
;Variable definitions
; These variables are only needed if low priority interrupts are used. 
; More variables may be needed to store other special function registers used
; in the interrupt routines.

		UDATA

WREG_TEMP	RES	1	;variable in RAM for context saving 
STATUS_TEMP	RES	1	;variable in RAM for context saving
BSR_TEMP	RES	1	;variable in RAM for context saving

		UDATA_ACS

NUM_TO_MUL_05	    RES	    1
NUM_TO_MUL_025	    RES	    1
NUM_TO_MUL_0125	    RES	    1	
NUM_TO_MUL_075	    RES	    1
NUM_TO_MUL_075_COPY RES	    1
 
Num_HighInt	    RES	    1
Num		    RES	    1
NumA		    RES	    1
NumA_prev	    RES	    1
NumB		    RES	    1	
NumB_prev	    RES	    1	
NumB_prev2	    RES	    1
NumA_RES	    RES	    1
NumA_prev_RES	    RES	    1
NumB_prev_RES	    RES	    1	
NumB_prev2_RES	    RES	    1
TMP_V		    RES	    1			    
		    
	

;******************************************************************************
;EEPROM data
; Data to be programmed into the Data EEPROM is defined here

DATA_EEPROM	CODE	0xf00000

		DE	"Test Data",0,1,2,3,4,5

;******************************************************************************
;Reset vector
; This code will start executing when a reset occurs.

RESET_VECTOR	CODE	0x0000

		goto	Main		;go to start of main code

;******************************************************************************
;High priority interrupt vector
; This code will start executing when a high priority interrupt occurs or
; when any interrupt occurs if interrupt priorities are not enabled.

HI_INT_VECTOR	CODE	0x0008

		bra	HighInt		;go to high priority interrupt routine

;******************************************************************************
;Low priority interrupt vector
; This code will start executing when a low priority interrupt occurs.
; This code can be removed if low priority interrupts are not used.

LOW_INT_VECTOR	CODE	0x0018

		bra	LowInt		;go to low priority interrupt routine

;******************************************************************************
;High priority interrupt routine
; The high priority interrupt code is placed here.

		CODE

HighInt:
;	*** high priority interrupt code goes here ***
		BCF INTCON,INT0IF; Сброс флага прерывания
		BSF LATB,RB2; Включение светодиода на RB2
		MOVLW .256; Запись 256 в WREG
		MOVWF Num_HighInt; Запись из WREG в перемнную Num_HighInt
		LFSR FSR0,0x100; Задаем базовый адрес для массива А - 0x100
		
Mint1:		MOVLW 0xF0; Запись 0xF0 в WREG
		MOVWF TMR1H; Запись из WREG в перемнную TMR1H
		MOVLW 0xC8; Запись 0xC8 в WREG
		MOVWF TMR1L; Запись из WREG в перемнную TMR1L
		BCF PIR1,TMR1IF; Сброс флага таймера
		BSF T1CON,TMR1ON; Включение таймера
Mint2:		BTFSS PIR1,TMR1IF; Задержка ~3.9 мс
		    BRA Mint2
		BSF ADCON0,1; Запуск АЦП
Mint3:		BTFSS ADCON0,1; Проверка на завершение оцифровки
		    BRA Mint3
		MOVFF ADRESH,TMP_V; Запись значений в массив А
		MOVLW 0x81; Запись 0x81 в WREG
		ADDWF TMP_V,F; Операция сложения WREG и TMP_V. Результат записывается в TMP_V
		MOVFF TMP_V,INDF0; Перенос значения из TMP_V в память МК
		INCF FSR0L,F; Переход на следующий адрес массива А
		DECFSZ Num_HighInt,F; Уменьшение на 1 количество опираций, и проверка на нулевое значение
		    BRA Mint1
		BSF LATB,RB3; Включение светодиода на RB3
		RETFIE FAST; Возравщение из прерывания 

;******************************************************************************
;Low priority interrupt routine
; The low priority interrupt code is placed here.
; This code can be removed if low priority interrupts are not used.

LowInt:
		movff	STATUS,STATUS_TEMP	;save STATUS register
		movff	WREG,WREG_TEMP		;save working register
		movff	BSR,BSR_TEMP		;save BSR register

;	*** low priority interrupt code goes here ***


		movff	BSR_TEMP,BSR		;restore BSR register
		movff	WREG_TEMP,WREG		;restore working register
		movff	STATUS_TEMP,STATUS	;restore STATUS register
		retfie

;******************************************************************************
;Start of main program
; The main program code is placed here.



Main:
	CALL CLRF_ALL; Вызов подпрограммы: Сброс всех переменных 
	NOP
	CALL CONF; Вызов подпрограммы: Конфигурация

M1:	BTFSS LATB,RB3; Ожидание включения светодиода на RB3
	    BRA M1
	
	MOVLW .256; Запись 256 в WREG
	MOVWF Num; Запись из WREG в переменную Num
	LFSR FSR0,0x100; Задаем базовый адрес для массива А - 0x100
	LFSR FSR1,0x300; Задаем базовый адрес для массива В - 0x300

M2:	MOVFF NumA,NumA_prev; Перенос значения А(n) в А(n-1)
	MOVFF NumB_prev,NumB_prev2; Перенос значения B(n-1) в B(n-2)
	MOVFF NumB,NumB_prev; Перенос значения B(n) в B(n-1)
	CLRF NumB; Обнуление значения B(n)
	
	MOVFF INDF0,NumA; Перенос значения из памяти МК в А(n)
	CALL MAIN_FUNC; Вызов подпрограммы: Обработка массива А
	MOVFF NumB,INDF1; Перенос значения B(n) в память МК
	
	INCF FSR1L,F; Переход на следующий адрес массива В
	INCF FSR0L,F; Переход на следующий адрес массива А
	DECFSZ Num,F; Уменьшение на 1 количество опираций, и проверка на нулевое значение
	    BRA M2
	
	BSF LATB,RB4; Включение светодиода на RB4
	
M3:	BTFSC PORTB,RB1; Ожидание нажатия кнопки на RB1
	    BRA M3
	   
	GOTO Main; Перезапуск основной программы 
  
CLRF_ALL:
	CLRF Num_HighInt; Сброс переменной Num_HighInt
	CLRF Num; Сброс переменной Num	
	CLRF NumA; Сброс переменной NumA
	CLRF NumA_RES; Сброс переменной NumA_RES
	CLRF NumA_prev; Сброс переменной NumA_prev
	CLRF NumA_prev_RES; Сброс переменной NumA_prev_RES
	CLRF NumB; Сброс переменной NumB
	CLRF NumB_prev; Сброс переменной NumB_prev
	CLRF NumB_prev_RES; Сброс переменной NumB_prev_RES
	CLRF NumB_prev2; Сброс переменной NumB_prev2
	CLRF NumB_prev2_RES; Сброс переменной NumB_prev2_RES
	CLRF NUM_TO_MUL_05; Сброс переменной NUM_TO_MUL_05
	CLRF NUM_TO_MUL_025; Сброс переменной NUM_TO_MUL_025
	CLRF NUM_TO_MUL_0125; Сброс переменной NUM_TO_MUL_0125
	CLRF NUM_TO_MUL_075; Сброс переменной NUM_TO_MUL_075
	CLRF NUM_TO_MUL_075_COPY; Сброс переменной NUM_TO_MUL_075_COPY
	CLRF TMP_V; Сброс переменной TMP_V
	RETURN; Возвращение в основную программу
	
CONF:
	;OSCILLATOR 
	MOVLW b'01101010'; Fosc = 4 MHz
	MOVWF OSCCON
	MOVLW b'10000000'; F = 8 MHz
	MOVWF OSCTUNE
	;ADC
	MOVLW b'00010100'; Включаем АЦП
	MOVWF ADCON2
	MOVLW b'00001110'; Аналоговый вход AN0, остальные - цифровые
	MOVWF ADCON1
	MOVLW b'00000001'; Используем канал AD0
	MOVWF ADCON0
	;TIMER
	MOVLW b'10000001';16-битный, Т=1мкс
	MOVWF T1CON
	
	BSF INTCON,INT0IE; Нажатие кнопки вызывает прерывания
	BCF INTCON2,INTEDG0; По отрицательному фронту
	BSF RCON,IPEN; Разрешение прерываний
	BSF INTCON,GIEH; Разрешение глобальных прерываний
	; Порты входных индикаторов
	BSF TRISA,RA0
	BSF TRISB,RB0
	BSF TRISB,RB1
	; Порты выходных индикаторов
	BCF TRISB,RB2
	BCF TRISB,RB3
	BCF TRISB,RB4
	BCF TRISB,RB5
	CLRF LATB; Сброс всех значений в регистре LATB
	
	RETURN; Возвращение в основную программу

		
MAIN_FUNC:
	MOVFF NumA,NUM_TO_MUL_0125; Запись значения NumA в переменную NUM_TO_MUL_0125
	CALL FUNC_MUL_TO_0125; Вызов подпрограммы: Умножение на 0.125
	MOVFF NUM_TO_MUL_0125,NumA_RES;
	
	MOVFF NumA_prev,NUM_TO_MUL_025; Запись значения NumA_prev в переменную NUM_TO_MUL_025
	CALL FUNC_MUL_TO_025; Вызов подпрограммы: Умножение на 0.25
	MOVFF NUM_TO_MUL_025,NumA_prev_RES; Запись значения NUM_TO_MUL_025 в переменную NumA_prev_RES
	
	MOVFF NumB_prev,NUM_TO_MUL_0125; Запись значения NumB_prev в переменную NUM_TO_MUL_0125
	CALL FUNC_MUL_TO_0125; Вызов подпрограммы: Умножение на 0.125
	MOVFF NUM_TO_MUL_0125,NumB_prev_RES; Запись значения NUM_TO_MUL_0125 в переменную NumB_prev_RES
	
	MOVFF NumB_prev2,NUM_TO_MUL_075; Запись значения NumB_prev2 в переменную NUM_TO_MUL_075
	CALL FUNC_MUL_TO_075; Вызов подпрограммы: Умножение на 0.75
	MOVFF NUM_TO_MUL_075,NumB_prev2_RES; Запись значения NUM_TO_MUL_075 в переменную NumB_prev2_RES
	
	MOVF NumA_RES,W; Запись значения NumA_RES в WREG
	ADDWF NumB,F; Операция сложения WREG и NumB. Результат записывается в NumB
	
	NEGF NumA_prev_RES; Перевод NumA_prev_RES в дополнительный код
	MOVF NumA_prev_RES,W; Запись значения NumA_prev_RES в WREG
	ADDWF NumB,F; Операция сложения WREG и NumB. Результат записывается в NumB
	
	MOVF NumB_prev_RES,W; Запись значения NumB_prev_RES в WREG
	ADDWF NumB,F; Операция сложения WREG и NumB. Результат записывается в NumB
	
	NEGF NumB_prev2_RES; Перевод NumA_prev_RES в дополнительный код
	MOVF NumB_prev2_RES,W; Запись значения NumB_prev2_RES в WREG
	ADDWF NumB,F; Операция сложения WREG и NumB. Результат записывается в NumB
	
	RETURN; Возвращение в основную программу
	
FUNC_MUL_TO_05:
    	BTFSS NUM_TO_MUL_05,7
	    BRA K1_05
	BRA K2_05
K1_05:	RRCF NUM_TO_MUL_05,F
	BCF NUM_TO_MUL_05,7
	BRA K3_05
K2_05:	RRCF NUM_TO_MUL_05,F
	BSF NUM_TO_MUL_05,7
	BRA K3_05
K3_05:	RETURN
	
FUNC_MUL_TO_025:
	BTFSS NUM_TO_MUL_025,7; Проверка старшего бита в переменной NUM_TO_MUL_025
	    BRA K1_025
	BRA K2_025
K1_025:	RRCF NUM_TO_MUL_025,F; Сдвиг переменной NUM_TO_MUL_025 на 1 разряд
	BCF NUM_TO_MUL_025,7; Обнуление старшего разряда
	RRCF NUM_TO_MUL_025,F; Сдвиг переменной NUM_TO_MUL_025 на 1 разряд
	BCF NUM_TO_MUL_025,7; Обнуление старшего разряда
	BRA K3_025
K2_025:	RRCF NUM_TO_MUL_025,F; Сдвиг переменной NUM_TO_MUL_025 на 1 разряд
	BSF NUM_TO_MUL_025,7; Установка старшего разряда
	RRCF NUM_TO_MUL_025,F; Сдвиг переменной NUM_TO_MUL_025 на 1 разряд
	BSF NUM_TO_MUL_025,7; Установка старшего разряда
	BRA K3_025
K3_025:	RETURN; Возвращение в основную программу
	
FUNC_MUL_TO_0125:
	BTFSS NUM_TO_MUL_0125,7; Проверка старшего бита в переменной NUM_TO_MUL_0125
	    BRA K1_0125
	BRA K2_0125
K1_0125:RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BCF NUM_TO_MUL_0125,7; Обнуление старшего разряда
	RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BCF NUM_TO_MUL_0125,7; Обнуление старшего разряда
	RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BCF NUM_TO_MUL_0125,7; Обнуление старшего разряда
	BRA K3_0125
K2_0125:RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BSF NUM_TO_MUL_0125,7; Установка старшего разряда
	RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BSF NUM_TO_MUL_0125,7; Установка старшего разряда
	RRCF NUM_TO_MUL_0125,F; Сдвиг переменной NUM_TO_MUL_0125 на 1 разряд
	BSF NUM_TO_MUL_0125,7; Установка старшего разряда
	BRA K3_0125
K3_0125:RETURN; Возвращение в основную программу
	
FUNC_MUL_TO_075:
	MOVFF NUM_TO_MUL_075,NUM_TO_MUL_025; Запись значения NUM_TO_MUL_075 в переменную NUM_TO_MUL_025
	CALL FUNC_MUL_TO_025; Вызов подпрограммы: Умножение на 0.25
	MOVFF NUM_TO_MUL_025,NUM_TO_MUL_075_COPY; Запись значения NUM_TO_MUL_025 в переменную NUM_TO_MUL_075_COPY
	NEGF NUM_TO_MUL_075_COPY; Перевод NUM_TO_MUL_075_COPY в дополнительный код
	MOVF NUM_TO_MUL_075_COPY,W; Запись значения NUM_TO_MUL_075_COPY в WREG
	ADDWF NUM_TO_MUL_075,F; Операция сложения WREG и NUM_TO_MUL_075. Результат записывается в NUM_TO_MUL_075
	RETURN; Возвращение в основную программу
	
	
;******************************************************************************
;End of program
		
    END


